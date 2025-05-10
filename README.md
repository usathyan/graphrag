# GraphRAG Querying Experiments: Detailed Report and Commentary

## 1. Introduction

This report details a series of experiments conducted to evaluate the capabilities of a GraphRAG instance in indexing a corpus of scientific articles (converted from PDF to Markdown) and answering complex domain-specific queries. The primary objectives were to:
*   Successfully build a knowledge graph from the processed articles.
*   Test the system's ability to retrieve and synthesize information for a range of scientific questions.
*   Observe the impact of different Large Language Models (LLMs) used for the chat/synthesis components on query success and answer quality.
*   Document the process, challenges, and outcomes to inform future use.

The experiments involved setting up the GraphRAG environment, processing input documents, indexing them, and then systematically running a predefined set of 10 queries using different configurations, primarily focusing on the `global_search` method.

## 2. Methodology

### 2.1. Environment Setup:
*   A Python virtual environment (`.venv`) was created using `uv`.
*   Dependencies were installed from `requirements.txt`.
*   The project was initialized as a Git repository.
*   Standard GraphRAG configuration files (`settings.yaml`, `prompts/`, etc.) were used and adapted.

### 2.2. Input Data:
*   The primary input data consisted of scientific articles, originally in PDF format, located in the `input/articles/` directory.
*   A utility script (`utils/pdf_to_markdown.py`) using the `markitdown` library was created and used to convert these PDFs into Markdown files, which were stored in `input/articles/processed/`.
*   The GraphRAG pipeline was configured to read these Markdown files.

### 2.3. Indexing Pipeline:
*   The GraphRAG indexing pipeline was run using the command: `graphrag index --root .` (after activating the virtual environment).
*   This process involved:
    *   Resetting previous outputs by deleting `output/` and `cache/` directories.
    *   Configuring `settings.yaml` to point `input:base_dir` to `input/articles/processed/` and `input:file_pattern` to `.*\\.md$$` to correctly identify the Markdown files.
    *   The indexing pipeline extracts entities and relationships, builds communities, generates summaries, and creates vector embeddings (stored in LanceDB by default).
*   The `default_embedding_model` was configured to use a local Ollama instance serving `nomic-embed-text`.

### 2.4. Querying Process:
*   A set of 10 predefined scientific queries was used.
*   Queries were executed using the command: `graphrag query --root . --method global --query "YOUR_QUERY_HERE" | cat`.
*   The `global_search` method was primarily used, which typically involves:
    1.  Mapping relevant document chunks to the query.
    2.  Reducing these mapped chunks into a synthesized answer, often involving LLM calls.

### 2.5. LLM Configurations for Chat/Synthesis (`default_chat_model` in `settings.yaml`):
The following models were used via OpenRouter API (`api_base: https://openrouter.ai/api/v1`), with the API key sourced from the `OPENROUTER_API_KEY` environment variable:
1.  `openai/gpt-4o`: Initial model. Encountered API credit exhaustion issues early in the querying phase.
2.  `openai/gpt-4.1-nano`: Switched to this model to potentially mitigate credit issues and test a smaller model.
3.  `openai/gpt-4o-mini`: Final model used for the most comprehensive set of query runs after adding more credits. This model offered a good balance of capability and (assumed) cost.

## 3. Results and Detailed Commentary

The majority of successful and detailed answers were obtained using the `openai/gpt-4o-mini` model after OpenRouter credits were replenished.

### 3.1. Performance of `openai/gpt-4o-mini` (Final Model):

*   **Query 1: "What are emerging therapeutic targets for non-small cell lung cancer identified in the last five years?"**
    *   **Result:** SUCCESSFUL.
    *   **Commentary:** Provided a comprehensive list including EML4-ALK, EGFR mutations, PTK7, and BRAF. It also correctly mentioned the role of machine learning and precision medicine in identifying these targets. The answer was well-structured and referenced data from the processed documents.

*   **Query 2: "Which proteins have been implicated as druggable targets in CRISPR screens for metabolic diseases?"**
    *   **Result:** FAILED ("I am sorry but I am unable to answer this question given the provided data.").
    *   **Commentary:** This query consistently failed across different models. It's highly probable that the indexed articles do not contain specific information about CRISPR screens for *metabolic diseases* identifying druggable protein targets.

*   **Query 3: "List novel kinase targets associated with resistance to current melanoma therapies."**
    *   **Result:** PARTIALLY SUCCESSFUL / NUANCED.
    *   **Commentary:** The model explained the context of resistance to BRAF/MEK inhibitors in melanoma. It stated that while kinase pathways are crucial, the provided documents did not *explicitly list novel kinase targets* for resistance. This indicates an ability to understand the query's nuance rather than just failing.

*   **Query 4: "Summarize recent advances in computational methods for target identification in rare genetic disorders."**
    *   **Result:** SUCCESSFUL.
    *   **Commentary:** Generated a detailed summary covering machine learning (GNNs, GCNs), multi-omics data integration, advanced models (deep autoencoders, TriModel, DeepWalk), and the role of collaborative efforts and databases (TTD, PharmGKB).

*   **Query 5: "What are the most frequently validated targets in published high-throughput screening studies for neurodegenerative diseases?"**
    *   **Result:** INCONSISTENT.
        *   First run with `gpt-4o-mini`: SUCCESSFUL. It identified voltage-gated ion channels (VGICs) and ligand-gated ion channels (LGICs).
        *   Second run (rerun of queries 1-7): FAILED ("I am sorry but I am unable to answer this question given the provided data.").
    *   **Commentary:** This inconsistency is important. It suggests that for queries where information might be sparsely represented or require significant synthesis, the LLM's path to an answer can vary between runs. The "map responses have score 0" warning in the failed run indicates the initial document mapping stage didn't find strong signals.

*   **Query 6: "Which disease pathways have newly identified protein targets with available structural data?"**
    *   **Result:** SUCCESSFUL.
    *   **Commentary:** Identified targets in cancer pathways (melanoma, TNBC, renal cell carcinoma, breast cancer, lung cancer) and cardiovascular disease pathways, mentioning specific genes/drugs and the role of databases like TTD/PDB.

*   **Query 7: "Find articles reporting on target deconvolution methods in phenotypic drug discovery."**
    *   **Result:** SUCCESSFUL.
    *   **Commentary:** Correctly identified the importance of these methods, referenced G.C. Terstappen's work in Nature Reviews Drug Discovery, and mentioned the role of ChEMBL and other supporting concepts.

*   **Query 8: "Summarize the use of knowledge graphs for predicting novel drug-target interactions."** (Run once with `gpt-4o-mini` after model switch)
    *   **Result:** SUCCESSFUL.
    *   **Commentary:** Provided a strong summary covering how KGs enhance predictive accuracy, integrate with machine learning (GNNs, GCNs), aid comprehensive understanding (off-target effects), and support drug repurposing.

*   **Query 9: "What are the most cited targets for immuno-oncology drug development in the last three years?"** (Run once with `gpt-4o-mini`)
    *   **Result:** FAILED ("I am sorry but I am unable to answer this question given the provided data.").
    *   **Commentary:** The specificity regarding "most cited" and "last three years" likely requires data (citation metrics, precise dating) not present or easily extractable from the narrative text of the articles.

*   **Query 10: "Which targets have been identified using multi-omics integration in cardiovascular disease research?"** (Run once with `gpt-4o-mini`)
    *   **Result:** FAILED ("I am sorry but I am unable to answer this question given the provided data.").
    *   **Commentary:** Similar to Query 9, the combination of "multi-omics integration" and "cardiovascular disease research" to yield a list of *specific targets* might not have been explicitly present in the indexed documents.

### 3.2. Observations on Model Comparison and Issues:

*   **`openai/gpt-4o` (Initial):** Performance could not be fully assessed due to rapid depletion of OpenRouter API credits. This highlights a critical operational aspect: ensure sufficient LLM credits for any substantial work.
*   **`openai/gpt-4.1-nano`:** This model was tested after the credit issue with `gpt-4o`. While it successfullyanswered some broader queries (e.g., Query 1, Query 4), it struggled with more specific ones, often returning "I am sorry but I am unable to answer..." or providing less detailed answers than `gpt-4o-mini`. This is expected given its smaller size and capability.
*   **`openai/gpt-4o-mini`:** This model provided the best balance of performance and (assumed) cost-efficiency for these experiments once credits were available. It successfully handled most of the complex summarization and information extraction queries.
*   **OpenRouter Credit Issue:** A recurring theme was the `openai.APIStatusError: Error code: 402 - {'error': {'message': 'Insufficient credits...'}}`. This appeared when the configured LLM on OpenRouter ran out of funds. It's crucial to monitor this when using paid API services. The errors appeared during the "map" phase of the global search, as this involves LLM calls to process and score document chunks.

### 3.3. General Challenges and Observations:
*   **Query Specificity:** Highly specific queries, especially those asking for ranked lists ("most cited") or quantitative details not typically found in narrative text, are challenging unless the input data is specifically structured or contains explicit mentions.
*   **Information Availability:** GraphRAG can only synthesize answers based on the information present in the indexed documents. If the information isn't there, no amount of prompt engineering or model power can create it.
*   **Inconsistency (Query 5):** The different outcomes for Query 5 on separate runs with the same model and data indicate the probabilistic nature of LLMs and perhaps variability in the document chunk retrieval and mapping steps. This is a factor to consider for applications requiring high determinism.
*   **Configuration (`settings.yaml`):** Small misconfigurations (e.g., incorrect `file_pattern` initially, needing `$$` to escape `$` for regex) can lead to pipeline failures. Careful configuration is key.

## 4. Conclusion and Recommendations

*   **Overall Success:** The GraphRAG system, particularly when paired with a capable model like `openai/gpt-4o-mini`, demonstrated a strong ability to index a corpus of scientific articles and provide detailed, synthesized answers to complex domain-specific questions.
*   **Model Choice:** `openai/gpt-4o-mini` appears to be a good candidate for this type of workload, offering a balance of capability and efficiency. Smaller models like `gpt-4.1-nano` may be suitable for simpler queries or when cost is an extreme constraint but will likely struggle with nuanced synthesis.
*   **Credit Management:** Active monitoring of API credits for services like OpenRouter is essential.
*   **Addressing Failed Queries:**
    *   For queries that failed due to likely absence of information (e.g., Query 2, 9, 10), the solution lies in curating a more appropriate or comprehensive set of input documents.
    *   For inconsistencies (e.g., Query 5), further investigation could involve analyzing the intermediate outputs of GraphRAG (if possible) or experimenting with different prompting strategies within the `prompts/` directory, or even different search methods (e.g., `local_search` if high precision on specific terms is needed, though `global_search` is generally better for broad questions).
*   **Future Work:**
    *   Experiment with different prompting strategies in `prompts/*.txt` for the global search steps (map, reduce) to see if answer quality or success on difficult queries can be improved.
    *   Explore the use of different embedding models if retrieval quality is a concern.
    *   For production systems, implement robust error handling and logging around API calls.

This set of experiments provides valuable insights into the practical application of GraphRAG and highlights areas for attention in configuration, model selection, and expectation management.

---
## Original Project Guide: Setup, Basic Usage, and Utilities

This section contains the original setup and usage instructions for this GraphRAG project example, including how to run the initial "A Christmas Carol" demonstration and use the PDF conversion utility.

For more in-depth information about GraphRAG itself, please also refer to the official resources:

👉 Official GraphRAG Repo: [microsoft/graphrag](https://github.com/microsoft/graphrag)  
👉 [Read the official docs](https://microsoft.github.io/graphrag/)  
👉 [GraphRAG Arxiv Paper](https://arxiv.org/abs/2404.16446)

### Overview (Original Project Example)

This repository provides a practical example of how to:
1. Set up a Python environment for GraphRAG.
2. Configure GraphRAG to use various LLM providers (including OpenAI, Azure OpenAI, OpenRouter.ai, and local models via Ollama).
3. Ingest local text data (`input/book.txt` in this example).
4. Run the GraphRAG indexing pipeline to build a knowledge graph, generate summaries, and create embeddings.
5. Query the indexed data using GraphRAG's global search capabilities.

### Setup Instructions

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
          type: openai # Ollama's OpenAI-compatible endpoint
          api_key: "ollama" # Or any non-empty string if Ollama doesn't require a key
          api_base: "http://localhost:11434/v1" # Default Ollama API endpoint
          model: "nomic-embed-text" # Ensure this model is pulled and served by Ollama
        ```
    *   Make sure model names match what's available on the service/Ollama instance.

### Running the Project (Original "A Christmas Carol" Example)

This example uses the text of "A Christmas Carol" (assumed to be in `input/book.txt`).

#### 1. Prepare Input Data

*   Place your desired text file(s) in the `input/` directory. For the demonstration, `input/book.txt` was used.

#### 2. Run the Indexing Pipeline

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

#### 3. Query the Index

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

### Important Considerations

*   **LLM Costs & Usage:** Be mindful of the costs associated with using proprietary LLMs (OpenAI, Azure OpenAI, OpenRouter). Indexing large datasets can consume significant tokens.
*   **Local Models (Ollama):** Ensure your Ollama instance is running and the specified models are downloaded (`ollama pull <model_name>`) before running GraphRAG.
*   **Configuration is Key:** The behavior and performance of GraphRAG heavily depend on the configurations in `settings.yaml` and the capabilities of the chosen LLMs.

### Prompt Tuning

For optimal results with your specific data, you may need to tune the prompts used by GraphRAG. The default prompts are located in the `prompts/` directory. Refer to the [official GraphRAG Prompt Tuning Guide](https://microsoft.github.io/graphrag/posts/prompt_tuning/overview/) for more details.

### Project Versioning

This project uses GraphRAG version `2.2.1` as specified in `requirements.txt`. For information on breaking changes or versioning of the core GraphRAG library, please refer to the official GraphRAG repository and its documentation.

### Utility Scripts

This project includes utility scripts located in the `utils/` directory.

#### PDF to Markdown Converter (`utils/pdf_to_markdown.py`)

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

**Example (as requested by user for scientific articles):**

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