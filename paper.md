---
title: 'pura: A Python package for cleaning chemical data quickly'
tags:
  - Python
  - chemistry
  - reactions
  - balancing
  - resolver
authors:
  - name: Kobi C. Felton
    orcid: 0000-0002-3616-4766
    affiliation: 1 
  - name: Adarsh Arun
    affiliation: "2, 1"
  - name: Alexei A. Lapkin
    orcid: 0000-0001-7621-0889
    corresponding: true # (This is how to denote the corresponding author)
    affiliation: "1, 2, 3" # (Multiple affiliations must be quoted)
affiliations:
 - name: Department of Chemical Engineering and Biotechnology, University of Cambridge, United Kingdom
   index: 1
 - name: Cambridge Centre for Advanced Research and Education in Singapore Ltd., Singapore
   index: 2
 - name: Innovation Centre in Digital
Molecular Technologies, Yusuf Hamied Department of Chemistry, University of Cambridge, United Kingdom
   index: 3

date: 1 March 2023
bibliography: paper.bib

---

# Summary

The last few years has a seen a boom in the use of machine learning in the field of organic chemistry. This includes models to predict the outcomes of single reactions, predict optimal reaction pathways and even predict the efficiency particular reaction recipes. These machine learning models often are trained on data collected from the academic and patent literature, which was not initially intended for cheminformatics purposes. 


# Statement of need

`pura` is a python package designed to clean messy chemical data quickly. It fills a gap in the existing open-source cheminformatics software landscape. There are now tools to automatically extract data from  unstrcutured sources[https://github.com/jiangfeng1124/ChemRxnExtractor, http://chemdataextractor.org/] and databases to store cleaned chemical data [https://github.com/open-reaction-database]. However, for many projects, the majority of time is spent on cleaning extracted data to prepare it for downstream tasks. 

Other things to cite: all the different models being built
https://scholar.google.com/citations?view_op=view_citation&hl=en&user=l015S80AAAAJ&citation_for_view=l015S80AAAAJ:7PzlFSSx8tAC

We identified to key tasks that many projects need: identifier resolution/standardization and reaction balancing. The former uses a combination of scraped database and asyncrhonous API calls to various servers to resolve names of chemicals to machine readdable formats. The later helps with finding the correct stoichiometric coefficients for reactions, an important step before building reaciton models.


# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

# Acknowledgements

We acknowledge contributions from Zhen Guo.

# References