## Project Overview
Xena users want to be able to run Gene Set Enrichment Analysis (GSEA) on the subgroups they created in Xena. To enable this functionality for Xena users, we built upon Xena's previous success with Differential Gene Expression analysis using the Appyter framework developed by the Ma’ayan lab [1]. In this project, we created a new Appyter tool that runs blitzGSEA [2] to offer Xena users GSEA. Previously, blitzGSEA was evaluated and found to have better performance than the traditional GSEA implementation while producing similar, biologically consistent results. 

## Appyter Background Information
Appyter is a software framework developed by researchers at the Ma’ayan lab. It allows users to run bioinformatic methods without having to interact with code directly. Users upload or fetch data and configure tweakable parameters on a web-based form; the variables inputted by the user on the form are inserted into a Jupyter Notebook, and the Notebook is executed. A Jupyter Notebook with the results is rendered for the user. The Notebook consists of markdown text and interactive/static plots along with the options to download the Notebook (in .ipynb, .html, or .zip format), view the source code, or run the Appyter locally. 

## Plotly & Dash Background Information
Plotly is an interactive library used to produce dozens of chart types, including 3D graphs, statistical and scientific charts, and SVG maps. The graphs can be viewed in Jupyter Notebooks, standalone HTML files, or integrated into Dash applications [5]. Dash is a low-code framework used for building web applications in Python, R, Julia, and F# [4]. It can tie interactive UI elements (dropdowns, sliders, and graphs) to Plotly Graphs, and developers can customize the sizing, positioning, colors, and fonts of elements to their liking. Both Plotly and Dash are open-source and MIT Licensed. 

## Objective
As many researchers and scientists utilize the Xena browser, consistency, accessibility, and functionality are extremely important; thus, the implementation of blitzGSEA Appyter on the Xena browser must be developed with that in mind. To ensure this, the project consisted of:
* Understanding the Appyter framework and building a “toy” Appyter using the tutorial provided by the Ma’ayan lab 
* Understanding a real-world Appyter example, such as the “Bulk RNA-seq analysis pipeline,” and deploying it
* Deploying the Xena browser on localhost 
* Understanding the Xena Differential Gene Expression Appyter
* Developing the Xena blitzGSEA Appyter

 ## Results & Outcomes
Throughout the first couple of weeks, a better understanding of the Appyter framework was gained through the examination of the existing Appyter Bulk RNA-seq Analysis Pipeline and deploying the Xena Differential Gene Expression (DE) Appyter connected to the Xena Browser on localhost. This allowed for the development of the Xena blitzGSEA Appyter to begin. 

The web form for the blitzGSEA Appyter (as seen in Figure 1) was developed in reference to the [Xena Differential Gene Expression Appyter](https://github.com/ucscXena/Xena_DE_Analysis_Pipeline), and many of the parameters were kept except those involved with Enrichr and L1000 small molecule analysis. Sections for selecting gene set libraries and blitzGSEA parameters were also added. The code from the last project on [blitzGSEA](https://github.com/callylin/xena_blitzGSEA) was implemented in addition to the DE Appyter code with some modifications to the production of the full table, top table, running sum plot, and detailed output. 

![webform](/images/xena_appyter_webform.png)
**Figure 1. blitzGSEA Appyter Web Form** — *The web form where users are able to set parameters for running blitzGSEA.*

The [source code for blitzGSEA](https://github.com/MaayanLab/blitzgsea) allows for a full table, a top table, and a detailed output and running sum plot for one gene set to be made. However, the interactivity of each is limited, as they are outputted in a static .tsv/.csv and .png format. To allow for better visualization and functionalities, the computing of the top table, running sum plot, and detailed output table were rewritten to utilize Dash (modified code available in Code Availability), a Python framework that enables the creation of interactive data visualization web apps [4]. As seen in Figure 2, users can now sort the data in ascending and descending order for a given column and utilize pagination to view all the data in one place. Furthermore, a button was added to allow users to download a .csv file for the detailed output. 

![do_table](/images/xena_appyter_do_table.gif)
**Figure 2. Plotly Dash Detailed Output Table** — *An example detailed output table in the blitzGSEA Appyter that uses pagination and has sorting for each column.*

The top table portion (as seen in Figure 3) of the Appyter is accessed through a link in the Notebook to a dashboard produced with Plotly Dash. Clicking on a gene set name produces an interactive running sum plot and detailed output for the gene set. While the graphs are being plotted, a circular spinner is displayed to indicate to the user that the graphs will appear momentarily. Originally, the dashboard was intended to be in the Notebook so users did not have to open a separate tab; however, upon further investigation and reaching out to one of the developers of Appyter (Daniel Clarke), it was found that due to the design of Appyter, it isn’t possible for end-users to send information to the kernel, as it would affect execution. This is why the dashboard is opened in a separate tab and the full table in the Notebook isn’t interactive.

![dashboard](/images/xena_appyter_dashboard.gif)
**Figure 3. Plotly Dash Top Table, Running Plot, and Detailed Output Dashboard** — *An example Dash app dashboard with a top table, running sum plot, and detailed output table. Upon clicking a gene set name, the plot and detailed output are updated accordingly.*

Overall, the blitzGSEA Appyter is now functional and provides researchers and scientists with the ability to run GSEA on their data more efficiently. In the next couple of months, the Appyter will be deployed in Docker and testing will be done to the user interface to ensure that researchers and scientists have an efficient and good experience analyzing their data. 
 
## References
[1] Clarke, D., Jeon, M., Stein, D. J., Moiseyev, N., Kropiwnicki, E., Dai, C., Xie, Z., Wojciechowicz, M. L., Litz, S., Hom, J., Evangelista, J. E., Goldman, L., Zhang, S., Yoon, C., Ahamed, T., Bhuiyan, S., Cheng, M., Karam, J., Jagodnik, K. M., Shu, I., … Ma'ayan, A. (2021). Appyters: Turning Jupyter Notebooks into data-driven web apps. Patterns (New York, N.Y.), 2(3), 100213. https://doi.org/10.1016/j.patter.2021.100213 

[2] Lachmann, A., Xie, Z., & Ma'ayan, A. (2022). blitzGSEA: Efficient computation of Gene Set Enrichment Analysis through Gamma distribution approximation. Bioinformatics (Oxford, England), btac076. Advance online publication. https://doi.org/10.1093/bioinformatics/btac076 

[3] Lab, M. (n.d.). What is an Appyter? Appyters. Retrieved June 24, 2022, from https://appyters.maayanlab.cloud/#/what-is-an-appyter/ 

[4] Dash documentation & user guide. Plotly. (n.d.). Retrieved August 25, 2022, from https://dash.plotly.com/   

[5] Getting Started with Plotly in Python. Plotly. (n.d.). Retrieved August 25, 2022, from https://plotly.com/python/getting-started/ 
