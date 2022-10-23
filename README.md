## Project Overview
Xena users want to be able to run Gene Set Enrichment Analysis (GSEA) on the subgroups they created in Xena. To enable this functionality for Xena users, we built upon Xena's previous success with Differential Gene Expression analysis using the Appyter framework developed by the Ma’ayan lab [1]. In this project, we created a new Appyter tool that runs blitzGSEA [2] to offer Xena users GSEA. Previously, blitzGSEA was evaluated to be significantly more performant than the traditional GSEA implementation while producing similar results that are biologically consistent. 

## Appyter Background Information
Appyter is a software framework developed by researchers at the Ma’ayan lab. It allows users to run bioinformatic methods without having to interact with code directly. Users upload or fetch data and configure tweakable parameters on a web-based form; the variables inputted by the user on the form are inserted into a Jupyter Notebook, and the Notebook is executed. A Jupyter Notebook with the results is rendered for the user. The Notebook consists of markdown text and interactive/static plots along with the options to download the Notebook (in .ipynb, .html, or .zip format), view the source code, or run the Appyter locally. 

## Objective
As many researchers and scientists utilize the Xena browser, consistency, accessibility, and functionality are extremely important; thus, the implementation of blitzGSEA Appyter on the Xena browser must be developed with that in mind. To ensure this, over the last ten weeks, the project consisted of:
* Understanding the Appyter framework and building a “toy” Appyter using the tutorial provided by the Ma’ayan lab 
* Understanding a real-world Appyter example, such as the “Bulk RNA-seq analysis pipeline,” and deploying it
* Deploying the Xena browser on localhost 
* Understanding the Xena Differential Gene Expression Appyter
* Developing the Xena blitzGSEA Appyter

 ## Results & Outcomes
Throughout the first couple of weeks, a better understanding of the Appyter framework was gained through the examination of existing Appyters (such as the Bulk RNA-seq Analysis Pipeline) and deploying the Xena Differential Gene Expression (DE) Appyter connected to the Xena Browser on localhost. This allowed for the development of the Xena blitzGSEA Appyter (accessible in Code Availability) to begin. 

The web form for the blitzGSEA Appyter (as seen in Figure 1) was developed in reference to the Xena Differential Gene Expression Appyter (available under Code Availability), and many of the parameters were kept beside those involved with Enrichr and L1000 small molecule analysis. Additionally, sections for selecting gene set libraries and blitzGSEA parameters were added. The code from the last project on blitzGSEA (as seen in Code Availability) was implemented in addition to the DE Appyter code with some modifications to the production of the full table, top table, running sum plot, and detailed output. 
