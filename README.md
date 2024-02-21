[![INFORMS Journal on Computing Logo](https://INFORMSJoC.github.io/logos/INFORMS_Journal_on_Computing_Header.jpg)](https://pubsonline.informs.org/journal/ijoc)

# The continuous time-resource tradeoff scheduling problem with time windows

This archive is distributed in association with the [INFORMS Journal on
Computing](https://pubsonline.informs.org/journal/ijoc) under the [MIT License](LICENSE). 
This is a code repo for experiments in the paper "The continuous time-resource tradeoff scheduling problem with time windows" by Christian Artigues,
Emmanuel Hébrard, Alain Quilliot, Hélène Toussaint.


## Cite

To cite the contents of this repository, please cite both the paper and this repo, using their respective DOIs.

https://doi.org/10.1287/ijoc.2022.0142

https://doi.org/10.1287/ijoc.2022.0142.cd

Below is the BibTex for citing this snapshot of the repository.

```
@misc{art2024,
  author =        {Christian Artigues and Emmanuel H\'ebrard and Alain Quilliot and H\'el\`ene Toussaint},
  publisher =     {INFORMS Journal on Computing},
  title =         {{The continuous time-resource tradeoff scheduling problem with time windows}},
  year =          {2024},
  doi =           {10.1287/ijoc.2022.0142},
  url =           {https://github.com/INFORMSJoC/2022.0142},
  note =          {Available for download at https://github.com/INFORMSJoC/2022.0142}
}  
```

## Instances

Instances are available in the data directory

The way instances have been generated is explained section 7.1 of the article.

## Replicating

- To run the code, you will need to make sure that you have already installed Cplex.
- You will need CMake for compilation.

To compile: 
```
cd src
mkdir build
cd build
cmake ..
make
```

This will create 3 executables : 
 - flowMip for running the flow-based exact method (see article section 5.1)
 - flowHeur for running the flow-based heuristic method (see article section 5.2)
 - BB for running the exact Branch & Bound method (see article section 6)


To execute:

Examples:

```
./flowMip 15 10 2 => run the flow-based exact method on all files in the 'data' directory whose name starts with 'gen1_15_10_2'(10 instances)
./flowHeur 15 10 2 => run the flow-based heuristic on all files in the 'data' directory whose name starts with 'gen1_15_10_2'(10 instances)
./BB 15 10 2 => run the Branch & Bound on all files in the 'data' directory whose name starts with 'gen1_15_10_2'(10 instances)
```

The results are stored in files starting by "res_"


## Support

Please contact [Christian Artigues](christian.artigues@laas.fr) or [Hélène Toussaint](helene.toussaint@uca.fr) if you have any questions.
