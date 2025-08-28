**CPSE_SIF: Coupled Photosynthesis, Stomatal conductance and Energy balance model utilizing the MLR-SIF approach to estimate photosynthesis.**

This is a new model that simulates photosynthesis using MLR-SIF (Mechanistic Light Reaction model using SIF as input) coupled with stomatal conductance and energy balance components.  This is built upon the existing coupled model of photosynthesis, stomatal conductance, and transpiration framework developed by Kim & Lieth (2003) (Kim & Lieth, 2003; Yang et al., 2009). Photosynthesis is simulated using the MLR SIF model by Gu et al. (2019). The Ball, Woodrow and Berry (BWB) model is used for stomatal conductance (Ball et al., 1987). The energy balance model is used for estimating leaf temperature as a function of stomatal conductance (boundary layer conductance), and other environmental variables (Campbell & Norman, 1998; Kim & Lieth, 2003). 

The three main modules in this model are photosynthetic gas exchange (Photosynthesis module), solar geometry (Solar module), and radiative transfer within the canopy (RadTrans module) to determine sunlit and shaded fractions. The Solar module estimates solar geometry and radiation data (e.g., solar elevation, azimuth, PAR, NIR, direct/diffuse fractions), and this information is used in the RadTrans module, which then simulates the partition of incoming light and leaf area into sunlit and shaded components. The solar geometry and radiative transfer components are derived from approaches described in Campbell & Norman (1998) and de Pury & Farquhar (1997). The photosynthetic module supports both the Farquhar-von Caemmerer-Berry (for C3) and Ball-Berry (for C4) carbon assimilation as well as MLR-SIF models. This enables crop models to utilize SIF when data is available, while maintaining FvCB functionality. This makes it easy for users to switch between methods depending on the data available.   


**References:** 

Gu, L., Han, J., Wood, J. D., Chang, C. Y.-Y., & Sun, Y. (2019). Sun-induced Chl fluorescence and its importance for biophysical modeling of photosynthesis based on light reactions. New Phytologist, 223(3), 1179–1191. https://doi.org/10.1111/nph.15796

Kim, S.-H., Sicher, R. C., Bae, H., Gitz, D. C., Baker, J. T., Timlin, D. J., & Reddy, V. R. (2006). Canopy photosynthesis, evapotranspiration, leaf nitrogen, and transcription profiles of maize in response to CO2 enrichment. Global Change Biology, 12(3), 588–600. https://doi.org/10.1111/j.1365-2486.2006.01110.x

Farquhar, G. D., von Caemmerer, S., & Berry, J. A. (1980). A biochemical model of photosynthetic CO2 assimilation in leaves of C3 species. Planta, 149(1), 78–90. https://doi.org/10.1007/BF00386231

von Caemmerer, S., & Farquhar, G. (1981). Some relationships between the biochemistry of photosynthesis and the gas exchange of leaves. Planta, 153(4). https://doi.org/10.1007/BF00384257

Han, J., Chang, C. Y.-Y., Gu, L., Zhang, Y., Meeker, E. W., Magney, T. S., Walker, A. P., Wen, J., Kira, O., McNaull, S., & Sun, Y. (2022). The physiological basis for estimating photosynthesis from Chla fluorescence. New Phytologist, 234(4), 1206–1219. https://doi.org/10.1111/nph.18045

Ball, J. T., Woodrow, I. E., & Berry, J. A. (1987). A Model Predicting Stomatal Conductance and its Contribution to the Control of Photosynthesis under Different Environmental Conditions. In J. Biggins (Ed.), Progress in Photosynthesis Research (pp. 221–224). Springer Netherlands. https://doi.org/10.1007/978-94-017-0519-6_48

Campbell, G. S., & Norman, J. M. (1977). An Introduction to Environmental Biophysics (G. S. Campbell & J. M. Norman, Eds.). Springer New York. https://doi.org/10.1007/978-1-4612-1626-1_1

De Pury, D. G. G., & Farquhar, G. D. (1997). Simple scaling of photosynthesis from leaves to canopies without the errors of big-leaf models. Plant, Cell and Environment, 20(5), 537–557. https://doi.org/10.1111/j.1365-3040.1997.00094.x


