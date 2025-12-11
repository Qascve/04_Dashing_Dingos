# Assignment 3: Groupwork 2

**Group Members:** Anaga Ambady, Chi Ching Cheng, Yian Liu, Zhou Yang  
**Submission Date:** 11/12/2025

---

## Individual Contributions

### [Anaga Ambady]
**Role(s):** Documentation Contributor, README Writer

**Specific Contributions:**  
- **Documentation:** Created `groupwork2/README.md` to document the project structure and overview.  
- **Documentation:** Updated and formatted `CONTRIBUTIONS.md` for groupwork2.  
- **Other Contributions:** Participated in group meetings, helped with documentation and formatting.  

**Key Commits:**  
- [46e9c5d] (2025-12-10) – Added `groupwork2/README.md` with project overview and structure.  
- [c1c773d] (2025-12-10) – Updated `CONTRIBUTIONS.md` for groupwork2 contributions.  
- [de47ac5] (2025-12-10) – Formatted `CONTRIBUTIONS.md` for better readability.  

---

### [Chi Ching Cheng (CCC)]
**Role(s):** Code Developer (Regression Analysis), Documentation Contributor

**Specific Contributions:**  
- **Code/Functions Written:** Developed `PP_Regress_loc.R` script for location-specific regression analysis.  
- **Code Functionality:** Implemented linear regression analysis separated by Location field, outputting results to CSV file.  
- **Documentation:** Created and maintained documentation files for groupwork2.  
- **Testing/Debugging Work:** Tested the regression script to ensure correct functionality.  

**Key Commits:**  
- [5044b83] (2025-11-25) – Created `PP_Regress_loc.R` script and initial regression results.

---

### [Yian Liu]
**Role(s):** Code Developer (Autocorrelation Analysis)

**Specific Contributions:**  
- **Code/Functions Written:** Developed initial version of `TAutoCorr.R` script for autocorrelation analysis.  
- **Code Functionality:** Implemented permutation test to assess significance of correlation between successive years' temperatures.  
- **Data Processing:** Added temperature data file (`KeyWestAnnualMeanTemperature.RData`).  
- **Testing/Debugging Work:** Tested the script and generated initial results including histogram visualization.  

**Key Commits:**  
- [e335743] (2025-11-28) – Created `TAutoCorr.R` script with permutation testing implementation, added data file and initial results.  

---

### [Zhou Yang]
**Role(s):** Technical Lead, GitHub Repository Manager, Code Enhancer

**Specific Contributions:**  
- **Code Enhancement:** Improved `TAutoCorr.R` script robustness and added LaTeX compilation support.  
- **Documentation:** Created LaTeX report (`florida.tex`) and PDF document with methods, results, and interpretation.  
- **Code Fixes:** Fixed file structure issues and corrected errors in `PP_Regress_loc.R`.  
- **Repository Management:** Set up git structure, managed branches, performed code reviews and merges.  
- **Testing/Debugging Work:** Performed final code testing and integration checks.  

**Key Commits:**  
- [be6a088] (2025-12-03) – Added LaTeX report (`florida.tex`), PDF document, compilation script (`compileLaTeX.sh`), and improved `TAutoCorr.R` robustness.  
- [9f019d9] (2025-12-04) – Fixed file structure and corrected errors in `PP_Regress_loc.R`.  
- [f7badcc] (2025-12-03) – Added git keep files for results directory.  
- [6a7cb2e] (2025-12-03) – Added git ignore configuration.  
- Multiple PR merges: #11, #12, #13, #16, #17 – Reviewed and merged pull requests.  

---

## Collaboration Process

**Meeting Schedule:**  
- Met online regularly for group work coordination.  
- Daily communication via WhatsApp for discussion and support.  

**Work Distribution Strategy:**  
- Two main tasks were divided among team members:  
  - **Yian Liu** took responsibility for the autocorrelation analysis task (`TAutoCorr.R`).  
  - **CCC** took responsibility for the regression analysis task (`PP_Regress_loc.R`).  
  - **Zhou Yang** handled technical setup, code enhancement, LaTeX documentation, and final integration.  
  - **Anaga Ambady** focused on documentation and README creation.  

**Key Decisions:**  
- Used separate branches (`florida`, `my-feature-branch`) for each major task before merging into `main`.  
- Ensured each member could work independently and test code before integration.  
- Decided to use LaTeX for formal report presentation as required.  

**Challenges and Solutions:**  
- Initial file structure needed reorganization ,resolved by Zhou through systematic refactoring.  
- PP_Regress script had errors that were identified and fixed during code review.  

---

## Code Review Summary

**Code Review Process:**  
- Pull requests and peer reviews were used to ensure code quality.  

**Pull Requests:**  
- **PR #11:** [merge main to florida] – Reviewed by: Zhou Yang  
- **PR #12:** [All group work final regression] – Reviewed by: Zhou Yang  
- **PR #13:** [florida] – Reviewed by: Zhou Yang  
- **PR #16:** [Florida] – Reviewed by: Zhou Yang  
- **PR #17:** [florida] – Reviewed by: Zhou Yang  

---

## Testing and Quality Assurance

- Each member tested their own code before pushing to respective branches.  
- Yian Liu tested `TAutoCorr.R` and verified permutation test results.  
- CCC tested `PP_Regress_loc.R` and verified CSV output format.  
- Zhou Yang performed final testing and integration checks before merging.  
- Verification was done after merging individual branches into `main`.  

---

## Declaration

We declare that the above contributions are accurate and that all team members participated actively in this group work.

**Signatures:**  
- Anaga Ambady: 11/12/2025  
- Chi Ching Cheng (CCC): 11/12/2025  
- Yian Liu: 11/12/2025  
- Zhou Yang: 11/12/2025
