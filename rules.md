# Project Governance & AI Rules (`rules.md`)

This document defines the strict, non-negotiable rules for any Artificial Intelligence (agent, copilot, or model) generating code, designing interfaces, or modifying files in the **BioMarker-R** project.

---

## 1. Aesthetic & UI Anti-Patterns (STRICTLY PROHIBITED)

Any AI assisting on this project must NEVER generate, suggest, or implement any of the following cliché design trends:

- ❌ **NO purple to blue gradients** (no indigo-to-violet linear fills, no neon cyan/purple combinations).
- ❌ **NO gradient hero text** (no `bg-clip-text text-transparent bg-gradient-to-r`).
- ❌ **NO emojis in headings or UI elements** (no 🚀, 🧬, 🔬, 🔥, 💡, 📊 in titles, labels, or navigation).
- ❌ **NO "Inter font everywhere"** default typography. Use a crisp, intentional pairing (e.g., editorial grotesque sans or clean monospaced data fonts like JetBrains Mono / IBM Plex Mono for figures and tables).
- ❌ **NO colored border cards** (no neon borders, no `border-cyan-500/50`).
- ❌ **NO glassmorphism cards** (no `backdrop-blur-md bg-white/10`, no frosted translucent sheets).
- ❌ **NO low-contrast dark mode** (no dark gray `#1a1a1a` on charcoal `#121212` with muted `#666` body text). Use high-contrast, accessible ratios ($> 7:1$).
- ❌ **NO 3 icon boxes in a row** (the cliché "Feature 1 - Feature 2 - Feature 3" marketing landing page pattern).
- ❌ **NO badge floating above the headline** (no pills stating "✨ v1.0 is now live" or "The Future of RNA-Seq").
- ❌ **NO Lucide icons everywhere** (do not plaster an icon next to every single noun or label).
- ❌ **NO untouched Shadcn UI** boilerplate defaults without bespoke styling and scientific customization.
- ❌ **NO fade-in on scroll** or intrusive intersection observer animations.
- ❌ **NO cursor-following beams**, glowing spotlights, or animated background canvas particle nets.
- ❌ **NO button opacity fades on hover** (use crisp solid state transitions or sharp border changes).
- ❌ **NO inconsistent spacing** (strictly stick to a cohesive 4px / 8px grid system: 4, 8, 12, 16, 24, 32, 48px).
- ❌ **NO em dashes everywhere** in copy or generated markdown.
- ❌ **NO generic passwords** or hardcoded mock credentials.
- ❌ **NO serif italic accents** inserted into sans-serif titles for artificial "classiness".
- ❌ **NO Space Grotesk + Instrument Serif** typography pairing (the overused "modern retro" tech aesthetic).
- ❌ **NO grain over a gradient** texture overlays.

---

## 2. Code Quality & Bioinformatics Pipeline Integrity Rules

1. **Deterministic Reproducibility**:
   - Always set seeds (`set.seed(42)`) before any stochastic step (e.g., PCA, permutation testing in GSEA).
   - Use `renv` to record exact package versions.
2. **Standard Bioconductor Workflow**:
   - Never run raw counts directly into PCA or clustering; always apply variance-stabilizing transformation (`vst()`) or regularized log transformation (`rlog()`).
   - Use raw integer count matrices for `DESeqDataSetFromMatrix()`. Never pass normalized counts (TPM, FPKM, or CPM) to `DESeq2`.
3. **Identifier Integrity**:
   - Always retain Ensembl gene IDs as permanent keys alongside HGNC gene symbols to avoid ambiguity from aliased or deprecated symbols.
4. **File Link Standards**:
   - In all documentation and responses, use GitHub-style markdown links with `file://` scheme when referencing local workspace files.

---

## 3. Mandatory Security & System Rules

1. **Zero Client Secrets**: Never include API keys, database credentials, or secret tokens in frontend or client JS.
2. **Zero Config Exposure**: Never commit `.env` or config files containing credentials to Git. Ensure `.gitignore` ignores all secret patterns.
3. **Clean Git History**: Never push secrets to git history. If an accidental commit occurs, rebase/purge immediately before pushing.
4. **Strict Database RLS**: When databases are integrated, never disable Row Level Security (RLS) or leave tables world-readable.
5. **Private Storage Buckets**: Never configure S3, Firebase, or Supabase buckets as public/world-writable.
6. **Rate Limiting Enforcement**: Mandatory rate limiting on all public API endpoints and authentication endpoints (e.g., via SlowAPI / Redis).
7. **Production Hardening**: Disable debug mode (`DEBUG = False`), OpenAPI docs in production if restricted, and developer tools in production builds.
8. **Protected Routes**: Never expose administrative or analytical compute-triggering routes without authentication and authorization checks.
9. **SQL & Query Sanitization**: Strictly prevent SQL, NoSQL, and command injection via parameter binding and Pydantic schema validation.
10. **Secure Password Hashing**: Passwords must never be stored in plain text or with weak algorithms (MD5/SHA1); utilize Argon2id or bcrypt with high work factors.
11. **Obfuscate Internal Errors**: Never show raw stack traces, file paths, or database error messages to end users in API responses.
12. **Server-Side Auth Authority**: Never rely exclusively on client-side authentication checks; every protected resource must be validated on the server.

---

## 4. The Three Highest-Risk Fixes First

To guarantee the security, stability, and integrity of BioMarker-R from day one, prioritize these three critical fixes above all others:

### 1. Fix 1: Eliminate Secrets from Frontend Bundles and Git Repositories
- **The Risk**: Committing `.env` files or exposing API keys (such as cloud database tokens, AI API keys, or private backend keys) in client-side Next.js code leaks credentials permanently to the public and scraping bots.
- **The Immediate Fix**:
  - Enforce a bulletproof `.gitignore` covering `.env`, `.env.*`, `*.pem`, `*.key`, `data/raw/`.
  - In Next.js, strictly forbid any prefix other than `NEXT_PUBLIC_` for variables intended for the browser.
  - The only client-side environment variable allowed is `NEXT_PUBLIC_API_BASE_URL`. All operational keys must remain securely on Render environment settings.

### 2. Fix 2: Server-Side Input Sanitization & Regex Boundary Validation
- **The Risk**: Malicious input in query parameters (e.g., `gene_symbol`, `search`, `fdr_cutoff`, `limit`) could cause directory traversal, denial of service (DoS) through unbounded memory consumption, or query injection.
- **The Immediate Fix**:
  - Implement strict Pydantic v2 schemas on every FastAPI route:
    - Gene symbols must match regex `^[A-Za-z0-9\-_]{1,30}$`.
    - Numeric parameters must enforce strict boundaries (e.g., $0.0 \le \text{FDR} \le 1.0$, $1 \le \text{limit} \le 1000$).
  - Reject unexpected query parameters and return standard HTTP 422/400 errors without internal debug tracebacks.

### 3. Fix 3: Strict Cross-Origin Resource Sharing (CORS) & Rate Limiting Enforcement
- **The Risk**: Allowing wildcard CORS (`allow_origins=["*"]`) combined with missing rate limiting leaves the Render backend vulnerable to unauthorized third-party embedding, data scraping, and compute exhaustion.
- **The Immediate Fix**:
  - Configure FastAPI `CORSMiddleware` with explicit origin whitelists:
    ```python
    allow_origins=[
        "https://biomarker-r.vercel.app",
        "http://localhost:3000"
    ]
    ```
  - Attach `SlowAPI` middleware to rate-limit all IP addresses (e.g., 60 requests/minute), shielding analytical endpoints from exhaustion attacks.
