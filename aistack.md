# Artificial Intelligence & LLM Stack (`aistack.md`)

## 1. Overview & Role of AI in BioMarker-R
In this project, Artificial Intelligence (AI) and Large Language Models (LLMs) are utilized strictly as **assistive reasoning, literature summarization, and development acceleration tools**. AI does not replace empirical statistical testing or biological validation.

### Core AI Tenet
> **Rule of Statistical Sovereignty**: AI models are strictly prohibited from generating, altering, imputing, or hallucinating raw count matrices, statistical test statistics ($t$-statistic, Wald test, quasi-likelihood), fold-change values, p-values, or survival hazard ratios. All numerical results originate exclusively from deterministic Bioconductor routines (`DESeq2`, `clusterProfiler`, `survival`).

---

## 2. AI Architecture & Integration Points

```
+-----------------------------------------------------------------------------------+
|                           AI ASSISTIVE CAPABILITIES                               |
|                                                                                   |
|  [Statistical Outputs] (DESeq2 / clusterProfiler)                                 |
|            |                                                                      |
|            v                                                                      |
|  [Deterministic Gene Shortlist] (Top DEGs: LFC > 2, FDR < 0.01, KM p < 0.05)      |
|            |                                                                      |
|            v                                                                      |
|  +-----------------------------------------------------------------------------+  |
|  |             Bio-Contextual LLM Synthesizer (Gemini / Claude)                |  |
|  |                                                                             |  |
|  |  Prompt Strategy: Grounded Literature Retrieval & Summarization             |  |
|  |  - Inputs: Gene Symbol (e.g., FOXA1), Log2FC (-3.4), Dysregulated Pathways   |  |
|  |  - References: NCBI Gene ID, UniProt Function, PubMed Primary Literature     |  |
|  |  - Task: Synthesize biological mechanism in TNBC oncogenesis                 |  |
|  |  - Output: Structured JSON conforming strictly to Pydantic BioSummary Schema|  |
|  +-----------------------------------------------------------------------------+  |
|            |                                                                      |
|            v                                                                      |
|  [Pre-Compiled Bio-Card Artifacts] (Cached JSON in backend, zero live LLM latency)|
+-----------------------------------------------------------------------------------+
```

---

## 3. Model Profiles & Use-Case Mapping

| Capability | Model Recommendation | Execution Mode | Scope & Deliverable |
| :--- | :--- | :--- | :--- |
| **Code Generation & Architecture** | Claude 3.5 Sonnet / Gemini 1.5 Pro | Pair Programming / IDE | Writing R Bioconductor scripts, FastAPI endpoints, Next.js components. |
| **Bioinformatics Domain Review** | Gemini 1.5 Pro (Deep Context) | Development Planning | Validating statistical model choices (e.g., negative binomial dispersion shrinkage, BH FDR adjustment). |
| **Gene Functional Summarization** | Gemini 1.5 Flash / Pro | Offline Script (Build Time) | Generates short, verified 2-paragraph summaries of top 20 candidate biomarkers in TNBC. |
| **Fast Literature Grounding** | PubMed E-Utilities + LLM RAG | Scripted / Offline | Links candidate biomarkers to published TNBC clinical trials and known mechanisms. |

---

## 4. Prompt Engineering & Grounding Standards

### Deterministic Output Enforcement
Whenever an LLM is used to synthesize biomarker biological context, the prompt must mandate:
1. **Explicit Identification**: Must cite approved HGNC symbol and Entrez ID.
2. **Hallucination Prevention**: If clinical association with TNBC is unconfirmed in literature, the LLM must explicitly state *"No direct TNBC mechanistic trial confirmed; association derived solely from transcriptomic survival stratification."*
3. **Structured Schema Validation**: Responses must validate against a strict Pydantic JSON schema:
   ```json
   {
     "gene_symbol": "STRING",
     "entrez_id": "INTEGER",
     "cellular_role": "STRING (max 50 words)",
     "tnbc_relevance": "STRING (max 80 words)",
     "key_pathways": ["ARRAY of STRING"],
     "literature_pmids": ["ARRAY of STRING"]
   }
   ```

---

## 5. Security & Deployment Safeguards for AI
1. **Zero Client-Side LLM Calls**: Web clients never call LLM provider APIs directly. No OpenAI or Google Gemini API keys are bundled into frontend JavaScript.
2. **Build-Time Generation**: All LLM biological summaries are pre-computed during data export (`06_export_artifacts.R` / helper Python script) and saved as static JSON. This eliminates runtime API billing, rate-limit failures, and latency in production.
3. **No Unsanitized User-to-LLM Pipelines**: The platform does not expose an open prompt box that forwards raw user inputs to an LLM, eliminating prompt injection and jailbreaking risks entirely.
