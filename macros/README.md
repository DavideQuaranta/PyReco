In this folder you will find the ROOT macros I used to perform the Fermi-Dirac fit on each strip. In the folder `data\timefit\`.

## Usage

To perform the time fit for a certain run, follow these steps:
- enter the `macros` folder and launch a ROOT shell by typing:
     ```sh
    root
- load the `apv_raw.c` macro:
    ```sh
    .L apv_raw.C
     ```
    this macro performs the fit and puts the results in a TTree. For each strip the fitted time and its error are saved.
- launch the `fit_apv_times.C` macro by typing:
    ```sh
    .x fit_apv_times.C
     ```
> Note: in the folder `data\timefit\` you will find an example file called `run2124_timefit.root` which was obtained by running this macro on the `run2124_xtalk1_xtalk2.root` file.