# Example request (perturbation design, `marker` mode)

We have scRNA-seq of mouse bone-marrow-derived macrophages, 3 wells treated with a cytokine for
24 h and 3 vehicle wells, about 20k cells, one library per well, integer counts in `X`, no
annotation. We want to know which genes respond in which macrophage sub-population, and whether
the agent finds the cytokine's canonical signalling read-outs without being told. Drugs are not
the question. Species is mouse.
