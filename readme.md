Note: being developed as we speak! 

# Extracting personalised latent dynamics using multilevel hidden Markov models
This webpage contains all the materials for an afternoon workshop on extracting personalised latent dynamics using multilevel hidden Markov models. The materials on this website are [CC-BY-4.0](https://creativecommons.org/licenses/by/4.0/) licensed.

![cc](https://mirrors.creativecommons.org/presskit/icons/cc.svg) ![by](https://mirrors.creativecommons.org/presskit/icons/by.svg)


## Course objectives

Facilitated by technological advances such as smartphones, smart watches, and sensors, it has become relatively easy and affordable to collect data on groups of individuals with a high temporal resolution: intensive longitudinal data (ILD). Due to the high sampling frequency, ILD can uniquely be used to study how psychological, behavioral and physiological processes unfold over time at the within-person level, and between person differences herein. When the dynamics over time of interest can be represented by a latent construct consisting of mutually exclusive categories, the hidden Markov model (HMM; Rabiner, 1989; de Haan-Rietdijk et al., 2017) is a promising novel approach. The HMM is a probabilistic, unsupervised, longitudinal machine learning method which uncovers empirically derived latent (i.e., hidden) states and the dynamics between these latent states over time. Utilising the multilevel framework, heterogeneity between individuals is accommodated, facilitating the study of individual specific dynamics and differences herein.

The workshop starts with a conceptual introduction on the (multilevel) hidden Markov model and how it fits together with ILD using an empirical example. This is followed by a hands-on workshop using the R CRAN package mHMMbayes. 

At the end of this session, participants have a firm grasp of the basics of the multilevel hidden Markov model, as well as the skills to start applying this method in their own work.


## Prerequisites

Please bring your laptops. We assume the following:

- you are comfortable with estimating and interpreting univariate and multivariate statistical models such as regression models
- you are familiar with the `R` programming language and you have a recent version installed
- it's a bonus if you are somewhat familiar with xx
- you have installed the following `R` packages on your computer:
  - `mHMMbayes` (version 1.0.0)
  - `ggplot2`

You can use the following code to install these at once:
```r
install.packages(c("ggplot2", "mHMMbayes"))
```
  

## Workshop schedule & materials

| Time  | Duration | Activity     | Content                                            | link |
| :---: | :------: | :----------- | :------------------------------------------------- | :--- |
| 13:30 | 45       | Lecture      | Introduction & multilevel hidden Markov model      | [`intro.pdf`](./lectures/01_introduction/intro.pdf) |
| 14:15 | 45       | Practical    | Fitting a mHMM + group level parameters            | [`intro.html`](./practicals/01_introduction/intro.html) |
| 15:00 | 15       | Break        |                                                    |      |
| 15:15 | 45       | Lecture      | Model selection and fit + subject level parameters | [`fit_and_subject_level.pdf`](./lectures/02_fit_and_subject_level/fit_and_subject_level.pdf) |
| 16:00 | 45       | Practical    | Model selection and fit + subject level parameters | [`fit_and_subject_level.html`](./practicals/02_fit_and_subject_level/fit_and_subject_level.html) |
| 16:45 | 15       | Conclusion   | Conclusion + questions                             |  [`conclusion.pdf`](./lectures/03_discussion/discussion.pdf)    |

You can download the dataset we have prepared from here: [`xx`](./data/xx.rds). Save it in a nicely accessible place, we will be using it in every practical.


## Additional links

- Tutorial vignette package: 
- Estimation vignette package: 
- HMM book zuchinni
- papers.. 


## Contact

For questions about this course, you can contact the instructor Emmeke ([e.aarts@uu.nl](mailto:e.aarts@uu.nl)) directly. 

## Development notes to me: 
- number of states (reference Phole 2017)
- obtaining and interpreting subject specific parameters 
    - checking state pattern similarity over subjects  
- model checking: 
    - convergence (rho and traceplots)
    - PPCs 

- discussion/final notes:  
    - how many observations do we need? 
    - please do sanity checks on your data! Visualize, and objective measures to detect shift in response (response times, no variance, etc)
    - multivariate data! 
    - covariates 


    

