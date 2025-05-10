# Latent Knowledge Graphs for Drug Discovery

## 1. Introduction

This report details a series of experiments conducted to evaluate the capabilities of a GraphRAG instance in indexing a corpus of scientific articles (converted from PDF to Markdown) and answering complex domain-specific queries. The primary objectives were to:
*   Successfully build a knowledge graph from the processed articles.
*   Test the system's ability to retrieve and synthesize information for a range of scientific questions.
*   Observe the impact of different Large Language Models (LLMs) used for the chat/synthesis components on query success and answer quality.
*   Document the process, challenges, and outcomes to inform future use.

The experiments involved setting up the GraphRAG environment, processing input documents, indexing them, and then systematically running a predefined set of 10 queries using different configurations, primarily focusing on the `global_search` method.

The experiments were conducted on a M1-Max Mac studio, running local models using ollama for indexing, and using OpenRouter.ai's models for querying.

## 2. Methodology

### 2.1. Environment Setup:
*   A Python virtual environment (`.venv`) was created using `uv`.
*   Dependencies were installed from `requirements.txt`.
*   The project was initialized as a Git repository.
*   Standard GraphRAG configuration files (`settings.yaml`, `prompts/`, etc.) were used and adapted.

### 2.2. Input Data:
*   The primary input data consisted of scientific articles, originally in PDF format, located in the `input/articles/` directory.
*   These documents were identified by Perplexity, along with the questions to ask against them. 
*   A utility script (`utils/pdf_to_markdown.py`) using the `markitdown` library was created and used to convert these PDFs into Markdown files, which were stored in `input/articles/processed/`.
*   The GraphRAG pipeline was configured to read these Markdown files.

### 2.3. Indexing Pipeline:
*   The GraphRAG indexing pipeline was run using the command: `graphrag index --root .` (after activating the virtual environment).
*   This process involved:
    *   Resetting previous outputs by deleting `output/` and `cache/` directories.
    *   Configuring `settings.yaml` to point `input:base_dir` to `input/articles/processed/` and `input:file_pattern` to `.*\\.md$$` to correctly identify the Markdown files.
    *   The indexing pipeline extracts entities and relationships, builds communities, generates summaries, and creates vector embeddings (stored in LanceDB by default).
*   The `default_embedding_model` was configured to use a local Ollama instance serving `nomic-embed-text`. The preference to use ollama was a cost based decision.

### 2.4. Querying Process:
*   A set of 10 predefined scientific queries was used.
*   Queries were executed using the command: `graphrag query --root . --method global --query "YOUR_QUERY_HERE" | cat`.
*   The `global_search` method was primarily used, which typically involves:
    1.  Mapping relevant document chunks to the query.
    2.  Reducing these mapped chunks into a synthesized answer, often involving LLM calls.

### 2.5. LLM Configurations for Chat/Synthesis (`default_chat_model` in `settings.yaml`):
The following models were used via OpenRouter API (`api_base: https://openrouter.ai/api/v1`), with the API key sourced from the `OPENROUTER_API_KEY` environment variable:
1.  `openai/gpt-4o`: Initial model - most expensive of the 3 models used.
2.  `openai/gpt-4.1-nano`: Switched to this model to test differences to quality, while being cheapest.
3.  `openai/gpt-4o-mini`: Final model used for the most comprehensive set of query runs. This model offered a good balance of capability and (assumed) cost.
4.  Using ollama models for querying resulted in errors due to incompability of outputs expected by GraphRAG library (JSON structured output responses), and modifying the library codebase was not the objective.

## 3. Results and Detailed Commentary

The majority of successful and detailed answers were obtained using the `openai/gpt-4o-mini` model.

### 3.1. Performance of `openai/gpt-4o-mini` (Final Model):

*   **Query 1: "What are emerging therapeutic targets for non-small cell lung cancer identified in the last five years?"**
    *   **Result:** SUCCESSFUL.
    *   **Commentary:** Provided a comprehensive list including EML4-ALK, EGFR mutations, PTK7, and BRAF. It also correctly mentioned the role of machine learning and precision medicine in identifying these targets. The answer was well-structured and referenced data from the processed documents.
    *   **Sample Response (using openai/gpt-4o-mini):**
        ```text
        The documents discuss several emerging therapeutic targets for non-small cell lung cancer (NSCLC) identified in recent years. These include:

        1.  **EML4-ALK Fusion Gene**: This is a well-established target, and drugs like crizotinib and alectinib are used. The documents mention its relevance in the context of resistance mechanisms and ongoing research.
        2.  **EGFR Mutations**: EGFR tyrosine kinase inhibitors (TKIs) are standard treatment, but resistance is a problem. The documents touch upon mechanisms of resistance and the search for new strategies against EGFR-mutant NSCLC.
        3.  **PTK7 (Protein Tyrosine Kinase 7)**: Identified as a target in various cancers, including NSCLC. It's implicated in Wnt signaling pathways and tumor progression.
        4.  **BRAF Mutations**: While more common in melanoma, BRAF mutations (like V600E) also occur in NSCLC and are druggable targets.
        5.  **MTBP (MDM2-Binding Protein)**: Mentioned as a potential novel target in NSCLC, with its overexpression linked to poor prognosis.

        The documents emphasize the role of precision medicine, understanding resistance mechanisms, and leveraging computational approaches (like machine learning and knowledge graphs) to identify and validate new therapeutic targets for NSCLC. The focus is often on overcoming resistance to existing therapies and finding targets in specific patient subpopulations.
        ```

*   **Query 2: "Which proteins have been implicated as druggable targets in CRISPR screens for metabolic diseases?"**
    *   **Result:** FAILED ("I am sorry but I am unable to answer this question given the provided data.").
    *   **Commentary:** This query consistently failed across different models. It's highly probable that the indexed articles do not contain specific information about CRISPR screens for *metabolic diseases* identifying druggable protein targets.
    *   **Sample Response (using openai/gpt-4o-mini):**
        ```text
        I am sorry but I am unable to answer this question given the provided data.
        ```

*   **Query 3: "List novel kinase targets associated with resistance to current melanoma therapies."**
    *   **Result:** PARTIALLY SUCCESSFUL / NUANCED.
    *   **Commentary:** The model explained the context of resistance to BRAF/MEK inhibitors in melanoma. It stated that while kinase pathways are crucial, the provided documents did not *explicitly list novel kinase targets* for resistance. This indicates an ability to understand the query's nuance rather than just failing.
    *   **Sample Response (using openai/gpt-4o-mini):**
        ```text
        Based on the provided documents, while specific novel kinase targets for melanoma therapy resistance aren't explicitly listed as a direct answer, the context suggests the following:

        *   **MAPK Pathway Kinases (BRAF, MEK):** Resistance to inhibitors targeting BRAF (like vemurafenib, dabrafenib) and MEK (like trametinib, cobimetinib) is a major challenge. Novel strategies often involve understanding and overcoming resistance mechanisms related to these or downstream/parallel kinases in the MAPK pathway.
        *   **Receptor Tyrosine Kinases (RTKs):** Aberrant RTK signaling can contribute to resistance. The documents mention PTK7 as an RTK, though its specific role in melanoma resistance isn't detailed as a novel kinase target *for resistance*.
        *   **PI3K/AKT/mTOR Pathway Kinases:** This pathway is often implicated in bypass resistance mechanisms to MAPK inhibitors.

        The documents emphasize that resistance to targeted therapies in melanoma is complex and can involve reactivation of the MAPK pathway, activation of bypass signaling pathways (often involving other kinases), or epigenetic modifications. The search for novel kinase targets is focused on those that mediate these resistance mechanisms.

        Therefore, while not a direct list, the provided information points towards continued investigation of kinases within the MAPK pathway, RTKs, and the PI3K/AKT/mTOR pathway as areas where novel targets associated with resistance might be identified or re-targeted.
        ```

*   **Query 4: "Summarize recent advances in computational methods for target identification in rare genetic disorders."**
    *   **Result:** SUCCESSFUL.
    *   **Commentary:** Generated a detailed summary covering machine learning (GNNs, GCNs), multi-omics data integration, advanced models (deep autoencoders, TriModel, DeepWalk), and the role of collaborative efforts and databases (TTD, PharmGKB).
    *   **Sample Response (using openai/gpt-4o-mini):**
        ```text
        Recent advances in computational methods for target identification in rare genetic disorders, as suggested by the provided documents, can be summarized as follows:

        1.  **Machine Learning and Graph-Based Techniques:**
            *   **Graph Neural Networks (GNNs), Graph Convolutional Networks (GCNs):** These are increasingly used to model complex biological relationships (gene-disease, protein-protein interactions) from heterogeneous data sources, helping to predict novel targets and understand disease mechanisms.
            *   **Knowledge Graphs (KGs):** KGs integrate diverse biological data (genomics, proteomics, literature, clinical data) to represent and infer relationships. They are used for predicting drug-target interactions (DTIs), identifying disease-associated genes, and drug repurposing. Methods like DeepWalk, Node2vec, and LINE are mentioned for learning embeddings from KGs.
            *   **Convolutional Neural Networks (CNNs):** Applied for tasks like predicting binding affinity and analyzing biological sequences.

        2.  **Integration of Multi-Omics Data:** Computational methods focus on integrating data from genomics, transcriptomics, proteomics, and metabolomics to build comprehensive models of disease and identify potential therapeutic targets. This holistic approach helps in understanding the complex interplay of factors in rare diseases.

        3.  **Advanced Predictive Models:**
            *   **Deep Autoencoders and Variational Autoencoders:** Used for feature learning and dimensionality reduction from high-dimensional biological data.
            *   **Ensemble Learning and Transfer Learning:** Applied to improve predictive accuracy, especially when data for specific rare diseases is scarce.
            *   Models like **TriModel** and techniques like **inductive matrix completion** are mentioned for DTI prediction.

        4.  **Network Medicine Approaches:** Analyzing biological networks (e.g., protein-protein interaction networks, gene regulatory networks) to identify key nodes (potential targets) and modules associated with rare genetic disorders.

        5.  **Leveraging Public Databases and Resources:** Computational methods extensively use public databases like OMIM, Orphanet, TTD (Therapeutic Target Database), DrugBank, and PharmGKB to gather information and validate findings.

        6.  **Focus on Explainability and Interpretability:** Efforts are being made to develop computational models that are not only predictive but also provide insights into the underlying biological mechanisms, which is crucial for target validation in rare diseases.

        The overall trend is towards more sophisticated AI/ML models that can handle diverse, large-scale biological data to uncover novel insights for diagnosing and treating rare genetic disorders by identifying and validating new therapeutic targets.
        ```

*   **Query 5: "What are the most frequently validated targets in published high-throughput screening studies for neurodegenerative diseases?"**
    *   **Result:** INCONSISTENT.
        *   First run with `gpt-4o-mini`: SUCCESSFUL. It identified voltage-gated ion channels (VGICs) and ligand-gated ion channels (LGICs).
        *   Second run (rerun of queries 1-7): FAILED ("I am sorry but I am unable to answer this question given the provided data.").
    *   **Commentary:** This inconsistency is important. It suggests that for queries where information might be sparsely represented or require significant synthesis, the LLM's path to an answer can vary between runs. The "map responses have score 0" warning in the failed run indicates the initial document mapping stage didn't find strong signals.
    *   **Sample Response (successful run with openai/gpt-4o-mini):**
        ```text
        Based on the provided documents, frequently validated targets in published high-throughput screening (HTS) studies for neurodegenerative diseases include:

        1.  **Voltage-gated ion channels (VGICs)**: These are mentioned as important drug targets for neurological disorders. The document "Computational Approaches for Drug Repurposing" highlights them in the context of drug discovery.
        2.  **Ligand-gated ion channels (LGICs)**: Similar to VGICs, LGICs are also indicated as significant targets for neurological conditions and are likely subjects of HTS studies.

        The documents also discuss the broader context of drug discovery for neurodegenerative diseases, mentioning:
        *   The use of HTS in identifying compounds that modulate disease-specific pathways.
        *   The importance of target identification and validation, which can be facilitated by computational methods and knowledge graphs.
        *   Challenges in CNS drug discovery, such as the blood-brain barrier penetration.

        While specific names of *most frequently validated* targets beyond these channel types are not exhaustively listed across all documents as a direct answer to HTS studies, the emphasis on ion channels as targets for neurological disorders is a recurring theme.
        ```

*   **Query 6: "Which disease pathways have newly identified protein targets with available structural data?"**
    *   **Result:** SUCCESSFUL.
    *   **Commentary:** Identified targets in cancer pathways (melanoma, TNBC, renal cell carcinoma, breast cancer, lung cancer) and cardiovascular disease pathways, mentioning specific genes/drugs and the role of databases like TTD/PDB.
    *   **Sample Response (using openai/gpt-4o-mini):**
        ```text
        Based on the provided documents, several disease pathways have newly identified protein targets with available structural data, often highlighted in the context of drug discovery and repurposing:

        1.  **Cancer Pathways (General and Specific):**
            *   **Melanoma:** Targets in the MAPK pathway (e.g., BRAF, MEK) have structural data, and ongoing research explores resistance mechanisms and new targets.
            *   **Triple-Negative Breast Cancer (TNBC), Renal Cell Carcinoma, Non-Small Cell Lung Cancer (NSCLC):** The documents refer to identifying targets for these cancers, and structural data is crucial for structure-based drug design (SBDD). Examples include kinases and proteins involved in cell cycle regulation or apoptosis.
            *   **PTK7 (Protein Tyrosine Kinase 7):** Implicated in various cancers (including NSCLC, breast cancer) and Wnt signaling. Structural information is valuable for designing inhibitors.

        2.  **Infectious Disease Pathways:**
            *   **Tuberculosis (Mycobacterium tuberculosis):** Targets like InhA, MmpL3, and other essential enzymes often have structural data available (e.g., from PDB) to guide drug development.

        3.  **Cardiovascular Disease Pathways:**
            *   The documents mention the use of multi-omics data and network analysis to identify novel targets in cardiovascular diseases. While specific *newly identified* protein targets with structural data aren't explicitly enumerated as a list, the approach implies that once identified, structural elucidation would be a priority for drug design.

        4.  **Neurodegenerative Disease Pathways:**
            *   Targets like voltage-gated and ligand-gated ion channels are mentioned, and structural data for these channel families is often available.

        **General Themes Regarding Structural Data:**
        *   The **Protein Data Bank (PDB)** is a key resource for available structural data.
        *   **Structure-Based Drug Design (SBDD)** heavily relies on the availability of 3D structures of protein targets.
        *   **Computational methods**, including homology modeling, are used to predict structures when experimental data is lacking, thus indirectly making structural information available for more targets.
        *   **Knowledge graphs** often integrate structural information from databases like PDB to aid in drug discovery and target validation.

        The documents emphasize that identifying a target and then obtaining its structural information (either experimentally or computationally) is a common workflow in modern drug discovery to facilitate the design of novel therapeutics.
        ```

*   **Query 7: "Find articles reporting on target deconvolution methods in phenotypic drug discovery."**
    *   **Result:** SUCCESSFUL.
    *   **Commentary:** Correctly identified the importance of these methods, referenced G.C. Terstappen's work in Nature Reviews Drug Discovery, and mentioned the role of ChEMBL and other supporting concepts.
    *   **Sample Response (using openai/gpt-4o-mini):**
        ```text
        The provided documents discuss target deconvolution methods in phenotypic drug discovery in several contexts:

        1.  **Importance and Definition:** Target deconvolution is highlighted as a crucial step after identifying a hit compound from a phenotypic screen. The goal is to understand the mechanism of action by identifying the specific molecular target(s) through which the compound exerts its effects. This is essential for further drug development and optimization.

        2.  **Computational Approaches:**
            *   **Knowledge Graphs and Databases:** Resources like ChEMBL, DrugBank, and the Therapeutic Target Database (TTD) are used. These databases contain information on known drug-target interactions, chemical structures, and bioactivity data, which can be computationally queried or used in machine learning models to predict potential targets for a hit compound.
            *   **Machine Learning:** Various ML techniques are employed, including those that predict drug-target interactions based on chemical structure, bioactivity profiles, or other features. These can help narrow down the list of potential targets.
            *   **Network Analysis:** Analyzing biological networks (e.g., protein-protein interaction networks) can help identify pathways perturbed by the compound, thus suggesting potential targets.

        3.  **Experimental Methods (Implied Context):** While the documents focus more on computational aspects and the role of knowledge graphs, they operate in the context of drug discovery pipelines that would also involve experimental target deconvolution techniques (e.g., affinity chromatography, proteomics-based methods like CETSA, genetic methods like CRISPR screens), though these experimental details are not the primary focus of the computational discussions.

        4.  **Examples and Mentions:**
            *   The article "A perspective on the future of drug discovery: the new purpose of AI and knowledge graphs" by G.C. Terstappen explicitly mentions phenotypic screening and the subsequent need for target identification.
            *   Discussions on drug repurposing often involve identifying new targets for existing drugs, which is a form of target deconvolution.

        In summary, the articles underscore the shift towards integrating computational methods, particularly those leveraging large-scale data and knowledge graphs, to facilitate and accelerate the target deconvolution process in phenotypic drug discovery. This helps in moving from an observed phenotype to a molecular understanding of drug action.
        ```

*   **Query 8: "Summarize the use of knowledge graphs for predicting novel drug-target interactions."** (Run once with `gpt-4o-mini` after model switch)
    *   **Result:** SUCCESSFUL.
    *   **Commentary:** Provided a strong summary covering how KGs enhance predictive accuracy, integrate with machine learning (GNNs, GCNs), aid comprehensive understanding (off-target effects), and support drug repurposing.
    *   **Sample Response (using openai/gpt-4o-mini):**
        ```text
        The provided documents extensively discuss the use of knowledge graphs (KGs) for predicting novel drug-target interactions (DTIs). Here's a summary of how KGs are utilized for this purpose:

        1.  **Enhanced Predictive Accuracy:**
            *   KGs integrate diverse biological data (e.g., drug information, protein/target data, disease associations, genomic data, textual evidence from literature) into a structured network format.
            *   This rich, heterogeneous information allows machine learning models, especially graph-based models (like GNNs, GCNs), to learn complex patterns and relationships that are not apparent from individual data sources alone.
            *   Methods like **DeepWalk, Node2vec, LINE, and TransE** are used to learn low-dimensional embeddings of entities (drugs, targets) and relations in the KG. These embeddings capture semantic similarities and can be fed into predictive models.
            *   Specific models mentioned for DTI prediction using KGs include **TriModel** (using tensor factorization) and various GNN architectures.

        2.  **Integration with Machine Learning:**
            *   KGs serve as a powerful input for various ML algorithms. Graph-based learning directly operates on the graph structure, while embeddings derived from KGs can be used in more traditional ML models (e.g., SVMs, random forests, deep neural networks).
            *   They help in addressing data sparsity issues often encountered in DTI prediction.

        3.  **Comprehensive Understanding and Interpretability:**
            *   KGs can represent not just direct interactions but also indirect relationships and pathways connecting drugs and targets.
            *   This allows for the prediction of novel DTIs based on shared pathways, similar drug properties, or analogous target functions.
            *   Some KG-based methods offer better interpretability by providing paths or subgraphs that explain a predicted interaction, which is valuable for validation.

        4.  **Drug Repurposing and Off-Target Prediction:**
            *   By identifying new targets for existing drugs (drug repurposing) or predicting potential off-targets, KGs contribute significantly to drug development and safety assessment.

        5.  **Incorporating Diverse Data Types:**
            *   KGs can naturally incorporate textual data (from biomedical literature via NLP), structural data (chemical structures of drugs, protein structures), omics data, and clinical trial information, providing a holistic view for DTI prediction.

        In essence, knowledge graphs provide a structured and semantically rich framework to represent and reason about complex biological systems. By applying advanced computational techniques on these KGs, researchers can more effectively predict novel drug-target interactions, leading to new therapeutic hypotheses and accelerating the drug discovery pipeline.
        ```

*   **Query 9: "What are the most cited targets for immuno-oncology drug development in the last three years?"** (Run once with `gpt-4o-mini`)
    *   **Result:** FAILED ("I am sorry but I am unable to answer this question given the provided data.").
    *   **Commentary:** The specificity regarding "most cited" and "last three years" likely requires data (citation metrics, precise dating) not present or easily extractable from the narrative text of the articles.
    *   **Sample Response (using openai/gpt-4o-mini):**
        ```text
        I am sorry but I am unable to answer this question given the provided data.
        ```

*   **Query 10: "Which targets have been identified using multi-omics integration in cardiovascular disease research?"** (Run once with `gpt-4o-mini`)
    *   **Result:** FAILED ("I am sorry but I am unable to answer this question given the provided data.").
    *   **Commentary:** Similar to Query 9, the combination of "multi-omics integration" and "cardiovascular disease research" to yield a list of *specific targets* might not have been explicitly present in the indexed documents.
    *   **Sample Response (using openai/gpt-4o-mini):**
        ```text
        I am sorry but I am unable to answer this question given the provided data.
        ```

### 3.2. Observations on Model Comparison and Issues:

*   **`openai/gpt-4o` (Initial):** Performance could not be fully assessed with this model during the initial phase.
*   **`openai/gpt-4.1-nano`:** This model was tested after the initial phase with `gpt-4o`. While it successfullyanswered some broader queries (e.g., Query 1, Query 4), it struggled with more specific ones, often returning "I am sorry but I am unable to answer..." or providing less detailed answers than `gpt-4o-mini`. This is expected given its smaller size and capability.It did well with summarization tasks.
*   **`openai/gpt-4o-mini`:** This model provided the best balance of performance and (assumed) cost-efficiency for these experiments. It successfully handled most of the complex summarization and information extraction queries. This model is also rated #1 on openrouter.ai's Science leaderboard.
*   It's crucial to monitor API usage when using paid API services. The errors appeared during the "map" phase of the global search, as this involves LLM calls to process and score document chunks.

### 3.2. Factual Graph Statistics and Examples (from processing 6 scientific articles)

The following statistics and examples were generated by running the `utils/inspect_graph_artefacts.py` script on the `output/` directory. The script prints available columns and then specific examples for each artefact type.

**1. Entities (`entities.parquet`)**
*   **Available Columns:** `['id', 'human_readable_id', 'title', 'type', 'description', 'text_unit_ids', 'frequency', 'degree', 'x', 'y']`
*   **Total Entities Extracted:** 3224 (inferred from full data, not shown in example)
*   **Commentary:** The `title` column appears to hold the primary human-readable name, and `description` provides further context.
*   **Example Entities (Top 5):**
    ```text
    1. Title: EUROPE PMC FUNDERS GROUP
       Description: A group involved in funding research, associated with the publication of manuscripts in the field of drug discovery
    2. Title: RITA SANTOS
       Description: Researcher involved in the study of molecular drug targets, affiliated with European Molecular Biology Laboratory-European Bioinformatics Institute
    3. Title: OLEG URSU
       Description: Researcher involved in the study of molecular drug targets, affiliated with the University of New Mexico School of Medicine
    4. Title: ANNA GAULTON
       Description: Researcher involved in the study of molecular drug targets, affiliated with European Molecular Biology Laboratory-European Bioinformatics Institute
    5. Title: A. PATRÍCIA BENTO
       Description: Researcher involved in the study of molecular drug targets, affiliated with European Molecular Biology Laboratory-European Bioinformatics Institute
    ```

**2. Relationships (`relationships.parquet`)**
*   **Available Columns:** `['id', 'human_readable_id', 'source', 'target', 'description', 'weight', 'combined_degree', 'text_unit_ids']`
*   **Total Relationships Identified:** 2242 (inferred from full data, not shown in example)
*   **Commentary:** The `source` and `target` columns contain entity titles, `description` explains the relationship, and `weight` indicates its strength/prominence.
*   **Example Relationships (Source -> Description -> Target, Weight, top 5):**
    ```text
    1. Source: EUROPE PMC FUNDERS GROUP
       Target: NAT REV DRUG DISCOV
       Description: The manuscript was published in Nat Rev Drug Discov with support from the Europe PMC Funders Group
       Weight: 7.0
    2. Source: RITA SANTOS
       Target: EUROPEAN MOLECULAR BIOLOGY LABORATORY-EUROPEAN BIOINFORMATICS INSTITUTE (EMBL-EBI)
       Description: Rita Santos is affiliated with EMBL-EBI
       Weight: 8.0
    3. Source: OLEG URSU
       Target: UNIVERSITY OF NEW MEXICO SCHOOL OF MEDICINE
       Description: Oleg Ursu is affiliated with the University of New Mexico School of Medicine
       Weight: 8.0
    4. Source: ANNA GAULTON
       Target: EUROPEAN MOLECULAR BIOLOGY LABORATORY-EUROPEAN BIOINFORMATICS INSTITUTE (EMBL-EBI)
       Description: Anna Gaulton is affiliated with EMBL-EBI
       Weight: 8.0
    5. Source: A. PATRÍCIA BENTO
       Target: EUROPEAN MOLECULAR BIOLOGY LABORATORY-EUROPEAN BIOINFORMATICS INSTITUTE (EMBL-EBI)
       Description: A. Patrícia Bento is affiliated with EMBL-EBI
       Weight: 8.0
    ```

**3. Community Reports (`community_reports.parquet`)**
*   **Available Columns:** `['id', 'human_readable_id', 'community', 'level', 'parent', 'children', 'title', 'summary', 'full_content', 'rank', 'rating_explanation', 'findings', 'full_content_json', 'period', 'size']`
*   **Total Community Reports Generated:** 167 (inferred from full data, not shown in example)
*   **Commentary:** The `summary` column provides a concise textual description of the community's theme.
*   **Top 2 Community Report Summaries (from 'summary' column):**
    ```text
    Report 1:
    The community centers around TriModel, an advanced knowledge graph embedding model designed for
    predicting drug-target interactions. TriModel is evaluated against several state-of-the-art models like NRLMF, BLM-NII, and COSINE, and consistently demonstrates superior performance. The model's development and evaluation are supported by various datasets and funding bodies, highlighting its significance in drug discovery and bioinformatics.

    Report 2:
    This community is centered around the publication of research on target deconvolution strategies in drug discovery by G.C. Terstappen in the journal Nature Reviews Drug Discovery. The relationship between the journal and the researcher highlights the significance of the research in the field of drug discovery.
    ```

These examples give a concrete indication of the knowledge graph constructed. The corpus size (6 documents) is relatively small, so a larger dataset would likely yield a more extensive and diverse graph.

**Inspecting Detailed Data:**
To explore the full content of these Parquet files, you can use the `utils/inspect_graph_artefacts.py` script (see instructions in the "Utility Scripts" section below) or load them directly using Python with libraries like `pandas` and `pyarrow`.

### 3.3. General Challenges and Observations:
*   **Query Specificity:** Highly specific queries, especially those asking for ranked lists ("most cited") or quantitative details not typically found in narrative text, are challenging unless the input data is specifically structured or contains explicit mentions.
*   **Information Availability:** GraphRAG can only synthesize answers based on the information present in the indexed documents. If the information isn't there, no amount of prompt engineering or model power can create it.
*   **Inconsistency (Query 5):** The different outcomes for Query 5 on separate runs with the same model and data indicate the probabilistic nature of LLMs and perhaps variability in the document chunk retrieval and mapping steps. This is a factor to consider for applications requiring high determinism.
*   **Configuration (`settings.yaml`):** Small misconfigurations (e.g., incorrect `file_pattern` initially, needing `$$` to escape `$` for regex) can lead to pipeline failures. Careful configuration is key.

## 4. Conclusion and Recommendations

*   **Overall Success:** The GraphRAG system, particularly when paired with a capable model like `openai/gpt-4o-mini`, demonstrated a strong ability to index a corpus of scientific articles and provide detailed, synthesized answers to complex domain-specific questions.
*   **Model Choice:** `openai/gpt-4o-mini` appears to be a good candidate for this type of workload, offering a balance of capability and efficiency. Smaller models like `gpt-4.1-nano` may be suitable for simpler queries or when cost is an extreme constraint but will likely struggle with nuanced synthesis.
*   Active monitoring of API usage for third-party LLM services is essential.
*   **Addressing Failed Queries:**
    *   For queries that failed due to likely absence of information (e.g., Query 2, 9, 10), the solution lies in curating a more appropriate or comprehensive set of input documents.
    *   For inconsistencies (e.g., Query 5), further investigation could involve analyzing the intermediate outputs of GraphRAG (if possible) or experimenting with different prompting strategies within the `prompts/` directory, or even different search methods (e.g., `local_search` if high precision on specific terms is needed, though `global_search` is generally better for broad questions).
*   **Future Work:**
    *   Experiment with different prompting strategies in `prompts/*.txt` for the global search steps (map, reduce) to see if answer quality or success on difficult queries can be improved.
    *   Explore the use of different embedding models if retrieval quality is a concern.
    *   For production systems, implement robust error handling and logging around API calls.

This set of experiments provides valuable insights into the practical application of GraphRAG and highlights areas for attention in configuration, model selection, and expectation management.

---
## Setup, Basic Usage, and Utilities

This section contains setup and usage instructions, including how to run the initial "A Christmas Carol" demonstration and use the PDF conversion utility.

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
*   **Local Models (Ollama):** Ensure your Ollama instance is running and the specified models are downloaded (`