# GraphRAG Project Example

This project demonstrates the use of Microsoft's GraphRAG to build a knowledge graph from text data and query it. It provides a self-contained, runnable example of setting up and executing the core GraphRAG indexing and querying workflow.

For more in-depth information about GraphRAG itself, please refer to the official resources:

👉 Official GraphRAG Repo: [microsoft/graphrag](https://github.com/microsoft/graphrag)  
👉 [Read the official docs](https://microsoft.github.io/graphrag/)  
👉 [GraphRAG Arxiv Paper](https://arxiv.org/abs/2404.16446)

## Overview

This repository provides a practical example of how to:
1. Set up a Python environment for GraphRAG.
2. Configure GraphRAG to use various LLM providers (including OpenAI, Azure OpenAI, OpenRouter.ai, and local models via Ollama).
3. Ingest local text data (`input/book.txt` in this example).
4. Run the GraphRAG indexing pipeline to build a knowledge graph, generate summaries, and create embeddings.
5. Query the indexed data using GraphRAG's global search capabilities.

## Setup Instructions

1.  **Prerequisites:**
    *   Git
    *   Python 3.10 or newer (this project was set up and tested with Python 3.11).
    *   `uv` (Python packaging tool). If you don't have `uv`, install it: `pip install uv`.

2.  **Clone the Repository (if you haven't already):**
    ```bash
    git clone <repository_url> # Replace with the URL of this repository
    cd <repository_directory_name>
    ```
    If you are already in the project directory, you can skip this step.

3.  **Create and Activate Virtual Environment:**
    ```bash
    uv venv --python python3.11 .venv # Or your preferred Python 3.10+ version
    source .venv/bin/activate
    ```

4.  **Install Dependencies:**
    ```bash
    uv pip install -r requirements.txt
    ```

5.  **Configure LLM and API Keys:**
    *   Create a `.env` file in the project root by copying `template.env` (if provided) or creating it from scratch. This file is `.gitignore`'d and should contain your API keys and other sensitive settings. 
        **Example `.env` content:**
        ```env
        # For OpenAI
        # GRAPHRAG_LLM_API_KEY="your_openai_api_key_here"
        # GRAPHRAG_LLM_MODEL="gpt-4-turbo-preview"

        # For Azure OpenAI
        # GRAPHRAG_LLM_API_KEY="your_azure_api_key"
        # GRAPHRAG_AZURE_OPENAI_ENDPOINT="https_your_azure_endpoint.openai.azure.com/"
        # GRAPHRAG_AZURE_OPENAI_DEPLOYMENT_NAME="your_deployment_name"
        # GRAPHRAG_LLM_MODEL="gpt-4-turbo-preview" # Or your Azure deployment model name

        # For OpenRouter.ai
        OPENROUTER_API_KEY="your_openrouter_key"

        # For local Ollama (usually no key needed for API access, but set model in settings.yaml)
        # OLLAMA_MODEL="llama3"
        ```
    *   Review and update `settings.yaml`. This file controls the behavior of GraphRAG. Pay close attention to the `llm`, `embeddings_llm`, and `community_reports_llm` sections to define the Language Models for different tasks.
        *   The `api_key` can reference environment variables like `${OPENROUTER_API_KEY}`.
        *   The `type` should be `openai`, `azure_openai`, or `static_response` (for testing).
        *   For OpenAI-compatible endpoints like OpenRouter or Ollama, you'll set `type: openai` and specify the `api_base` and `model`.

6.  **Using with OpenRouter.ai or Ollama:**
    *   Ensure the relevant API key (for OpenRouter) is in your `.env` file.
    *   In `settings.yaml`, modify the LLM configuration blocks. 

        **Example for OpenRouter (main LLM tasks):**
        ```yaml
        llm:
          type: openai
          api_key: ${OPENROUTER_API_KEY}
          api_base: "https://openrouter.ai/api/v1"
          model: "mistralai/mistral-7b-instruct" # Or your desired OpenRouter model
          # Other parameters like max_tokens, temperature as needed
        ```

        **Example for local Ollama (embeddings):**
        ```yaml
        embeddings_llm:
          type: openai # Ollama\'s OpenAI-compatible endpoint
          api_key: "ollama" # Or any non-empty string if Ollama doesn't require a key
          api_base: "http://localhost:11434/v1" # Default Ollama API endpoint
          model: "nomic-embed-text" # Ensure this model is pulled and served by Ollama
        ```
    *   Make sure model names match what's available on the service/Ollama instance.

## Running the Project

This example uses the text of "A Christmas Carol" (assumed to be in `input/book.txt`).

### 1. Prepare Input Data

*   Place your desired text file(s) in the `input/` directory. For the demonstration, `input/book.txt` was used.

### 2. Run the Indexing Pipeline

*   Ensure your virtual environment is active: `source .venv/bin/activate`
*   From the project root, run:
    ```bash
    graphrag index --root .
    ```
*   This command processes data in `input/`, builds the knowledge graph, and stores artifacts in `output/` and `cache/`. This can be time-consuming and may incur LLM costs.

    *Expected output snippet during indexing:*
    ```
    🚀 LLM Config Params Validated
    🚀 Embedding LLM Config Params Validated
    Running standard indexing.
    🚀 create_base_text_units
    ...
    🚀 generate_text_embeddings
    ...
    🚀 All workflows completed successfully.
    ```

### 3. Query the Index

*   After indexing, query your data:
    ```bash
    graphrag query --root . --method <method_name> --query "<your_question>"
    ```
*   Supported methods: `global`, `local`, `drift`, `basic`.

*   **Example Global Query (as demonstrated):**
    ```bash
    graphrag query --root . --method global --query "What is this book about?" | cat
    ```

*   **Example Output for "A Christmas Carol":**
    ```
    SUCCESS: Global Search Response:
    ### Overview
    
    The book revolves around the character of Ebenezer Scrooge, a miserly and solitary old man whos
    e life undergoes a profound transformation following a series of ghostly visits...
    (...rest of the example output...)
    ```

## Important Considerations

*   **LLM Costs & Usage:** Be mindful of the costs associated with using proprietary LLMs (OpenAI, Azure OpenAI, OpenRouter). Indexing large datasets can consume significant tokens.
*   **Local Models (Ollama):** Ensure your Ollama instance is running and the specified models are downloaded (`ollama pull <model_name>`) before running GraphRAG.
*   **Configuration is Key:** The behavior and performance of GraphRAG heavily depend on the configurations in `settings.yaml` and the capabilities of the chosen LLMs.

## Prompt Tuning

For optimal results with your specific data, you may need to tune the prompts used by GraphRAG. The default prompts are located in the `prompts/` directory. Refer to the [official GraphRAG Prompt Tuning Guide](https://microsoft.github.io/graphrag/posts/prompt_tuning/overview/) for more details.

## Project Versioning

This project uses GraphRAG version `2.2.1` as specified in `requirements.txt`. For information on breaking changes or versioning of the core GraphRAG library, please refer to the official GraphRAG repository and its documentation.

## Utility Scripts

This project includes utility scripts located in the `utils/` directory.

### PDF to Markdown Converter (`utils/pdf_to_markdown.py`)

This script converts PDF files to Markdown format using the `markitdown` library. It can process a single PDF file or all PDF files within a specified directory.

**Prerequisites:**

Ensure `markitdown` and its dependencies for PDF processing are installed. If you followed the main setup and installed `requirements.txt`, this should already be handled. The relevant line in `requirements.txt` is:
`markitdown[pdf,docx,pptx,xlsx]==0.1.1`

**Usage:**

Ensure your virtual environment is active (`source .venv/bin/activate`).

*   **To convert a single PDF file:**
    ```bash
    python3 utils/pdf_to_markdown.py <path_to_input_pdf_file>
    ```
    The Markdown file will be saved in the same directory as the input PDF with a `.md` extension.
    To specify a different output directory for the single converted file:
    ```bash
    python3 utils/pdf_to_markdown.py <path_to_input_pdf_file> -o <path_to_output_directory>
    ```

*   **To convert all PDF files in a directory:**
    ```bash
    python3 utils/pdf_to_markdown.py <path_to_input_directory> -o <path_to_output_directory>
    ```
    The `<path_to_output_directory>` is required when processing a directory. Markdown files will be saved in this output directory.

**Example (as requested by user):**

To convert all PDF files in `input/articles/` and save the Markdown output to `input/articles/processed/`:

1.  Create the input directory and place some PDF files in it (if not already present):
    ```bash
    mkdir -p input/articles
    # Add your PDF files to input/articles/
    ```
2.  Run the script:
    ```bash
    python3 utils/pdf_to_markdown.py input/articles -o input/articles/processed
    ```

## Scientific Article Querying Experiments (Recent Example)

Beyond the initial "A Christmas Carol" example, a series of experiments were conducted to test GraphRAG's capabilities on a corpus of scientific articles. This section summarizes that process and its outcomes.

**1. Objective:**
To index a collection of scientific articles (converted from PDF to Markdown) and evaluate GraphRAG's ability to answer complex, domain-specific queries using different LLMs for response synthesis.

**2. Input Data:**
*   Scientific articles were first converted from PDF to Markdown using the `utils/pdf_to_markdown.py` script.
*   The processed Markdown files were stored in `input/articles/processed/`.
*   GraphRAG was configured to use this directory as the `input:base_dir` with `file_pattern: ".*\\.md$$"` in `settings.yaml`.

**3. Indexing:**
The standard indexing pipeline (`graphrag index --root .`) was run to process these articles. The `default_embedding_model` used was `nomic-embed-text` served via a local Ollama instance.

**4. Querying and LLM Models Tested:**
A set of 10 scientific queries were run using the `global` search method (`graphrag query --root . --method global --query "..."`).
The following chat/synthesis LLMs (configured as `llm:model` in `settings.yaml` via OpenRouter) were tested:
*   `openai/gpt-4o`: Initial tests were impacted by API credit exhaustion.
*   `openai/gpt-4.1-nano`: Showed some success but was generally less capable for complex synthesis.
*   `openai/gpt-4o-mini`: This model yielded the most comprehensive and successful results once API credits were sufficient.

**5. Key Outcomes & Observations:**
*   **Successful Queries:** With `openai/gpt-4o-mini`, many queries related to therapeutic targets, computational methods, and knowledge graph applications in drug discovery yielded detailed and relevant answers based on the indexed articles.
*   **Failed Queries:** Some queries, particularly those asking for highly specific data not likely present in the narrative text of the articles (e.g., "most cited targets," information on CRISPR screens for *metabolic diseases* specifically), returned an "unable to answer" response. This highlights that GraphRAG's success is contingent on the information actually being present in the source documents.
*   **Model Impact:** The choice of LLM for the synthesis step significantly impacted the quality and detail of the answers. `gpt-4o-mini` provided a good balance.
*   **API Credit Management:** The experiments underscored the critical importance of ensuring sufficient API credits when using third-party LLM providers like OpenRouter. Errors due to `Insufficient credits` were encountered.
*   **Configuration:** Correctly configuring `settings.yaml` (e.g., `input:file_pattern`, LLM endpoints, and model names) is crucial for the pipeline to run successfully.

**Refer to the detailed experimental report and discussion (available in the development/chat history where these experiments were conducted) for a query-by-query breakdown and further commentary.**

**Reminder on Configuration:**
*   Ensure your `OPENROUTER_API_KEY` (or other relevant keys) is set in your `.env` file.
*   Carefully configure the `llm`, `embeddings_llm`, etc., sections in `settings.yaml` to point to your desired models and API endpoints. 