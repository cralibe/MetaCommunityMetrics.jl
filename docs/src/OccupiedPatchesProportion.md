# Occupied Patches Proportion Function
```@meta
CurrentModule = MetaCommunityMetrics
```
The Occupied Patches Proportion function in `MetaCommunityMetrics` provides a simple yet powerful metric for understanding the distribution and prevalence of species across different habitat patches within a metacommunity. By calculating the proportion of patches occupied by each species, this function helps ecologists assess the spatial extent of species distributions and identify potential patterns of rarity or commonness across the landscape.

This function draws on the concepts discussed by Ehrlén & Eriksson (2000) in their study on dispersal limitation and patchy occupancy in forest herbs. According to their findings, low occupancy may indicate dispersal limitation or strong competition, while high occupancy could suggest mass effects due to high dispersal rates or the ability to thrive in various conditions.

## Functionality Overview
The Occupied Patches Proportion metric quantifies the proportion of habitat patches in which a species is present. This information is crucial for understanding species distributions, particularly in fragmented landscapes or patchy environments where species may not occupy all available habitat. The OPP metric can be used to identify widespread species, which occupy a large number of patches, as well as rare species, which are restricted to only a few patches.

After calculating the proportion of patches occupied for each species, the species’ mean, minimum, and maximum occupancy can be derived. These values are valuable indicators of ecological processes:

- Low Occupancy: May indicate dispersal limitation or strong competition among species. Such patterns could suggest that certain species struggle to colonize or persist in many patches.
- High Occupancy: Could point to mass effects, where species are abundant in many patches, possibly due to high dispersal rates or the ability to thrive across a range of conditions.

## The Function
```@docs
prop_patches
```
