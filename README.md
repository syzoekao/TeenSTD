# Using multiple outcomes of sexual behavior to provide insights into chlamydia transmission and the effectiveness of prevention interventions in adolescents

This repository contains an `R` package (`TeenSTD`) for the published study, "Using multiple outcomes of sexual behavior to provide insights into chlamydia transmission and the effectiveness of prevention interventions in adolescents." Mathematical models developed for evaluating interventions are usually developed for the disease of interest alone, and calibrated to the observed epidemiological trend (e.g., disease incidence) of the disease of interest. However, incorporating additional outcomes generated by the same underlying disease dynamics may improve the precision of the unknown parameters estimated in the calibration process and lead to different policy implication. In this study, we modeled the spread of chlamydia among heterosexual adolescents aged 15-19 years in Minnesota in two calibration scenarios: (1) **chlamydia calibration**: calibrating the model to chlamydia incidence alone from 2005 to 2013, and (2) **dual calibration**: calibration the model to both chlamydia incidence and pregnancy trend from 2005-2013. The details of the study objective, methods, and findings can be found [here](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5711443/pdf/nihms875380.pdf).

The model code maintained in this repository was developed for the **dual calibration** scenario, which calibrated to both the chlamydia incidence and pregnancy incidence in 2005-2013 in Minnesota. The model code and the related datasets were provided in the `R` package. In this `README` file, we summarized the structure of the repository, and the components and the use of the package. Here is the outline of this `README` document: 

1. Provided the instructions regarding downloading and installing the `TeenSTD` package.
2. Outlined the structure of the package (raw data, model parameters, and model code). 
3. Executed the model code using calibrated parameters. 

## Installing `TeenSTD` package

Before downloading the folder and installing the package, please make sure that `R` is installed on your computer or laptop. If you have `R` installed, the following steps showed how to download and install the packages. 

* Download/clone the repository


    * If you have terminal, open terminal and type: 
    
    ```
    git clone https://github.com/syzoekao/TeenSTD.git </code>
    ```
    
    * If you don't have terminal, please download the zip file 
















