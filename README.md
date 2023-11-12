# GCA
The Guiding Center Approximation (GCA) code utilizes standard relativistic equations as outlined in Northrop's work from 1963. To apply this method effectively, it's essential that the characteristic time and spatial scales of particle gyromotion are significantly smaller than those of magnetic field variations.

This code employs the Bulirsch-Stoer method to solve particle motion within pre-specified analytical electric and magnetic fields. The method is known for its high accuracy and efficiency, particularly in cases where the fields exhibit smooth variations.

In this implementation, particles operate based on the Guiding Center approximation method. The primary goal is to determine the distribution function within a 3D magnetic field, which includes oscillations and an external dipolar field loaded from a separate file. The code leverages parallel computing capabilities, making it well-suited for scenarios involving numerous test particles.

This particular version of the code is designed for continuous runs over different time intervals, typically 25 TimeSteps. Users can specify a collection box to gather information about both particles that remain within it and those passing through it. This feature is commonly utilized to analyze density variations within specific loops and boxes.

The GCA code is improved and optimized. In this update, I've transitioned from Fortran 77 to Fortran 90, benefitting from enhanced features and improved readability. The code structure has been maintained without altering the logic or functionality. The adoption of Fortran 90's array syntax, explicit loop control, and modular design contributes to a more modern and maintainable codebase. Additionally, the code now adheres to Fortran 90 standards, making it compatible with modern compilers and facilitating future modifications and advancements in scientific computing.

# Attribution
Researchers who use this GCA code for scientific research are asked to cite the following two papers by Mehdi Yousefzadeh listed below.

1. Mehdi Yousefzadeh, Yao Chen, Hao Ning, Mahboub Hosseinpour, “Harmonic Electron Cyclotron Maser Emission along the
Coronal Loop”, The Astrophysical Journal, IOP science. DOI: https://doi.org/10.3847/1538-4357/ac6de3

2. Mehdi Yousefzadeh, Hao Ning, Yao Chen, “Harmonic Electron Cyclotron Maser Emission Excited by Energetic Electrons
Traveling inside a Coronal Loop”, The Astrophysical Journal, IOP science. DOI: https://doi.org/10.3847/1538-4357/abd8d5
