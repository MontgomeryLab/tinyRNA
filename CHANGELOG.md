### 11/21/2021
- The following improvements have been made to the  class_charts plot:
    - The proportion calculation: 
      - was: `(normalized class count) / (total normalized assigned counts)`
      - now: `(raw class count) / (total mapped reads)`
    - Due to the change in calculation, an "Unassigned" category has been added
    - A table with percentage proportion per class has been added
- The pipeline now supports low-DF experiment designs. If the experiment's degrees of freedom is less than 1, DESeq2 will refuse analysis. This condition is detected at the beginning of the pipeline, and if met, the DGE step will be skipped.

### 11/11/2021

- The following Alignment Stats have been corrected. **NOTE**: These issues are only relevant within the context of the Alignment Stats file. The issues and their fixes have no affect on outputs from Counter that are used in downstream analyses.:
  - Assigned Single-Mapping Reads
  - Assigned Multi-Mapping Reads
  - Sequences Assigned to Multiple Features
  
- The outdated and unclear description of the Features Sheet has been corrected
- Resource link for WS279 added to TUTORIAL.md
- Conda installation files now reside in their own top-level directory
- Bugfix: DGE-class scatter plots now show `Class=unknown` groups by default
- Plot stylesheet changed to match title font size with axes labels
- CHANGELOG.md created 