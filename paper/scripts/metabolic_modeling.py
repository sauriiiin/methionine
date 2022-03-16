import cobra
from cobra import Reaction, Metabolite
import math

# Load model
model = cobra.io.read_sbml_model('/Users/zoydwheeler1/Desktop/Carvunis_Lab_2/methionine/metabolic_modeling_2/ecYeastGEM_batch.xml')

# Correct metabolite ids:
for met in model.metabolites:
    met.id = met.id.replace('__91__', '_')
    met.id = met.id.replace('__93__', '')

#Solve the default model
summary = model.optimize()
print(summary)
print()

#Remove existing YLL reaction(s)
model.remove_reactions(['r_0815No1'],True)
model.remove_reactions(['arm_r_0815'],True)
model.remove_reactions(['r_0815_REVNo1'],True)
# model.remove_reactions(['arm_r_0815_REV'],True)

#Solve the model
summary = model.optimize()
print(summary)
print()

#Add YLL reaction
reaction = Reaction('YLL_rxn')
reaction.name = 'YLL058W_hypothetical_reaction'
reaction.lower_bound = 0
reaction.upper_bound = math.inf
hydrogen_sulfide = model.metabolites.get_by_id('s_0841[c]')
O_acetyl_homoserine = model.metabolites.get_by_id('s_1233[c]')
homocysteine = model.metabolites.get_by_id('s_1012[c]')
acetate = model.metabolites.get_by_id('s_0362[c]')
hydrogen = model.metabolites.get_by_id('s_0794[c]')
YLL = model.metabolites.get_by_id('prot_Q12198[c]')
reaction.add_metabolites({
    hydrogen_sulfide: -1.0,
    O_acetyl_homoserine: -1.0,
    YLL: -2.31481E-06,
    homocysteine: 1.0,
    acetate: 1.0,
    hydrogen: 1.0

})
print(reaction.reaction)
reaction.gene_reaction_rule = 'YLL058W'
print(reaction.genes)
model.add_reactions([reaction])
print()

#Add YLL reverse reaction
reaction2 = Reaction('YLL_rxn_REV')
reaction2.name = 'YLL058W_hypothetical_reaction_REV'
reaction2.lower_bound = 0
reaction2.upper_bound = math.inf
reaction2.add_metabolites({
    hydrogen_sulfide: 1.0,
    O_acetyl_homoserine: 1.0,
    YLL: -2.31481E-06,
    homocysteine: -1.0,
    acetate: -1.0,
    hydrogen: -1.0
})
print(reaction2.reaction)
reaction2.gene_reaction_rule = 'YLL058W'
print(reaction2.genes)
model.add_reactions([reaction2])
print()

#Solve the model
summary = model.optimize()
print(summary)
print()

#Remove MET15
cobra.manipulation.remove_genes(model,['YLR303W'])
print()

#Solve the model
summary = model.optimize()
print(summary)
print()







