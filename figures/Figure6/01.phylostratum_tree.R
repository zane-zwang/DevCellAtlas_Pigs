# Purpose: phylostratum tree.

library(ggplot2)
library(dplyr)
library(ggsci)
library(scales)
library(ggtree)
library(ape)

ordered_phylostratum <- sc.data %>%
  arrange(rank) %>%
  pull(phylostratum) %>%
  unique()

tree_text = "((((((((((((((((((((((cellular_organisms,Eukaryota),Opisthokonta),Metazoa),Eumetazoa),Bilateria),Deuterostomia),Chordata),Vertebrata),Gnathostomata),Euteleostomi),Sarcopterygii),Dipnotetrapodomorpha),Tetrapoda),Amniota),Mammalia),Theria),Eutheria),Boreoeutheria),Laurasiatheria),Artiodactyla),Suidae),Sus_scrofa);"

tree <- read.tree(text = tree_text)

ggtree(tree, layout = "rectangular", branch.length = 'none') +
  geom_tiplab(aes(label = label), hjust = -0.1, size = 3) +
  geom_nodepoint(size=2, shape = 16, color = 'gray')+
  ggtitle("Evolutionary Tree with Gene Counts") +
  theme_tree2()
