# immune-response-to-dengue
Stochastic model of the immune response to dengue infections
This model is designed to simulate the human response of B cell and T cells to the dengue virus. 
The immune system is modeled as a system of chemical reactions using the Gillespie algorithm (Woo & Reifman, Proc Nat Acad Sci USA, 2012). 
It uses the immune shape space model developed by Smith et al. (Smith et al., J Theor Biol 1997) and the immune response model by B cells for malaria by Chaudhury et al. (Journal of Immunology. 2014). 

REQUIREMENTS:

The immune modeling code requires Python 2.4 or later.

USAGE EXAMPLES: 

python run_dengueinfections.py <data_file> monovalent
python run_dengueinfections.py <data_file> polyvalent

CITATION:

Nguyen et al. Stochastic models of the adaptive immune response predict disease severity and captures enhanced cross-reactivity in natural dengue infections. Journal of Immunology. 2021
