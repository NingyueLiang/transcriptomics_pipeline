#!/usr/bin/env python3
"""
HVD API wrapper for calling Harvard's OpenAI Direct API.

This module provides a simple interface for making requests to the Harvard
OpenAI Direct API endpoint.
"""

import logging
import os
import time
from typing import Optional

try:
    import requests
except ImportError:
    requests = None


def call_hvd_api(
    prompt: str,
    model: str = "gpt-5-mini",
    api_key: Optional[str] = None,
    api_url: Optional[str] = None,
    max_retries: int = 5,
    timeout: int = 120,
) -> str:
    """
    Call the Harvard OpenAI Direct API (HVD API) to generate a response.
    
    Args:
        prompt: The prompt text to send to the model.
        model: Model name to use (default: "gpt-5-mini").
        api_key: API key for authentication. If None, will try to read from env var HVD_API_KEY.
        api_url: API URL. If None, uses default Harvard endpoint.
        max_retries: Maximum number of retry attempts (default: 5).
        timeout: Request timeout in seconds (default: 120).
    
    Returns:
        The generated text response from the model.
    
    Raises:
        ImportError: If requests library is not installed.
        ValueError: If API key is not provided.
        RuntimeError: If API call fails after all retries.
        ValueError: If no text content is found in the response.
    """
    if requests is None:
        raise ImportError("requests is not installed. Please install it with: pip install requests")
    
    # Get API key from parameter or environment variable
    if api_key is None:
        api_key = os.environ.get("HVD_API_KEY")
        if api_key is None:
            raise ValueError(
                "HVD API key not provided. Please set HVD_API_KEY environment variable "
                "or pass api_key parameter."
            )
    
    # Set API URL
    if api_url is None:
        api_url = "https://go.apis.huit.harvard.edu/ais-openai-direct/v1/responses"
    
    # Prepare payload for HVD API
    payload = {
        "model": model,
        "input": [
            {
                "role": "user",
                "content": prompt,
            },
        ],
    }
    
    # Prepare headers
    headers = {
        "Content-Type": "application/json",
        "api-key": api_key,
    }
    
    # Retry logic with exponential backoff
    base_delay = 1.0  # seconds
    
    logging.info(f"Requesting response from Harvard API (model: {model})")
    
    for attempt in range(max_retries):
        try:
            response = requests.post(
                api_url,
                headers=headers,
                json=payload,
                timeout=timeout
            )
            
            # Check for HTTP errors
            if response.status_code == 429:  # Rate limit
                if attempt < max_retries - 1:
                    delay = base_delay * (2 ** attempt)
                    logging.warning(
                        f"Rate limited. Retrying after {delay:.1f}s "
                        f"(attempt {attempt + 1}/{max_retries})..."
                    )
                    time.sleep(delay)
                    continue
                else:
                    raise RuntimeError(f"Rate limited after {max_retries} attempts")
            elif response.status_code == 503:  # Service unavailable
                if attempt < max_retries - 1:
                    delay = base_delay * (2 ** attempt)
                    logging.warning(
                        f"Service unavailable. Retrying after {delay:.1f}s "
                        f"(attempt {attempt + 1}/{max_retries})..."
                    )
                    time.sleep(delay)
                    continue
                else:
                    raise RuntimeError(f"Service unavailable after {max_retries} attempts")
            
            response.raise_for_status()
            response_data = response.json()
            
            # Extract text from HVD response format
            # Response format: {"output": [{"type": "message", "content": [{"type": "output_text", "text": "..."}]}]}
            output = response_data.get("output", [])
            content = ""
            
            for item in output:
                if item.get("type") == "message" and "content" in item:
                    for content_item in item["content"]:
                        if content_item.get("type") == "output_text" and "text" in content_item:
                            content = content_item["text"]
                            break
                    if content:
                        break
            
            if not content:
                raise ValueError("No text content found in HVD API response")
            
            logging.info("Received API response (%d characters).", len(content))
            return content
            
        except requests.exceptions.HTTPError as e:
            # Print the error response for debugging
            try:
                error_detail = response.json()
                logging.error(f"HTTP Error {response.status_code}: {error_detail}")
            except:
                logging.error(f"HTTP Error {response.status_code}: {response.text}")
            
            if attempt < max_retries - 1:
                delay = base_delay * (2 ** attempt)
                logging.warning(
                    f"Request error: {e}. Retrying after {delay:.1f}s "
                    f"(attempt {attempt + 1}/{max_retries})..."
                )
                time.sleep(delay)
            else:
                raise RuntimeError(f"Request error after {max_retries} attempts: {e}")
        except requests.exceptions.Timeout:
            if attempt < max_retries - 1:
                delay = base_delay * (2 ** attempt)
                logging.warning(
                    f"Request timeout. Retrying after {delay:.1f}s "
                    f"(attempt {attempt + 1}/{max_retries})..."
                )
                time.sleep(delay)
            else:
                raise RuntimeError(f"Request timeout after {max_retries} attempts")
        except requests.exceptions.RequestException as e:
            if attempt < max_retries - 1:
                delay = base_delay * (2 ** attempt)
                logging.warning(
                    f"Request error: {e}. Retrying after {delay:.1f}s "
                    f"(attempt {attempt + 1}/{max_retries})..."
                )
                time.sleep(delay)
            else:
                raise RuntimeError(f"Request error after {max_retries} attempts: {e}")
        except Exception as e:
            logging.error(f"Unexpected error calling HVD API: {e}")
            raise

