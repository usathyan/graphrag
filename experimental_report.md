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
*   **`openai/gpt-4.1-nano`:** This model was tested after the credit issue with `gpt-4o`. While it successfully answered some broader queries (e.g., Query 1, Query 4), it struggled with more specific ones, often returning "I am sorry but I am unable to answer..." or providing less detailed answers than `gpt-4o-mini`. This is expected given its smaller size and capability.
*   **`openai/gpt-4o-mini`:** This model provided the best balance of performance and (assumed) cost-efficiency for these experiments once credits were available. It successfully handled most of the complex summarization and information extraction queries.
*   **OpenRouter Credit Issue:** A recurring theme was the `openai.APIStatusError: Error code: 402 - {'error': {'message': 'Insufficient credits...'}}`. This appeared when the configured LLM on OpenRouter ran out of funds. It's crucial to monitor this when using paid API services. The errors appeared during the "map" phase of the global search, as this involves LLM calls to process and score document chunks.

### 3.3. General Challenges and Observations:
*   **Query Specificity:** Highly specific queries, especially those asking for ranked lists ("most cited") or quantitative details not typically found in narrative text, are challenging unless the input data is specifically structured or contains explicit mentions.
*   **Information Availability:** GraphRAG can only synthesize answers based on the information present in the indexed documents. If the information isn't there, no amount of prompt engineering or model power can create it.
*   **Inconsistency (Query 5):** The different outcomes for Query 5 on separate runs with the same model and data indicate the probabilistic nature of LLMs and perhaps variability in the document chunk retrieval and mapping steps. This is a factor to consider for applications requiring high determinism.
*   **Configuration (`settings.yaml`):** Small misconfigurations (e.g., incorrect `file_pattern` initially, needing `$$` to escape `$` for regex) can lead to pipeline failures. Careful configuration is key.

## 4. Conclusion and Recommendations

*   **Overall Success:** The GraphRAG system, particularly when paired with a capable model like `openai/gpt-4o-mini`, demonstrated a strong ability to index a corpus of scientific articles and provide detailed, synthesized answers to complex domain-specific questions.
*   **Model Choice:** `gpt-4o-mini` appears to be a good candidate for this type of workload, offering a balance of capability and efficiency. Smaller models like `gpt-4.1-nano` may be suitable for simpler queries or when cost is an extreme constraint but will likely struggle with nuanced synthesis.
*   **Credit Management:** Active monitoring of API credits for services like OpenRouter is essential.
*   **Addressing Failed Queries:**
    *   For queries that failed due to likely absence of information (e.g., Query 2, 9, 10), the solution lies in curating a more appropriate or comprehensive set of input documents.
    *   For inconsistencies (e.g., Query 5), further investigation could involve analyzing the intermediate outputs of GraphRAG (if possible) or experimenting with different prompting strategies within the `prompts/` directory, or even different search methods (e.g., `local_search` if high precision on specific terms is needed, though `global_search` is generally better for broad questions).
*   **Future Work:**
    *   Experiment with different prompting strategies in `prompts/*.txt` for the global search steps (map, reduce) to see if answer quality or success on difficult queries can be improved.
    *   Explore the use of different embedding models if retrieval quality is a concern.
    *   For production systems, implement robust error handling and logging around API calls.

This set of experiments provides valuable insights into the practical application of GraphRAG and highlights areas for attention in configuration, model selection, and expectation management. 