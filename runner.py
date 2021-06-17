from rdkit.Chem import MolFromSmiles as convert
from utils import metrics
from prob import Problem
from evol import Evolution




with open('moses/dataset/data/dataset.smi') as f:
    content = f.readlines()
# you may also want to remove whitespace characters like `\n` at the end of each line
smiles = [x.strip() for x in content] 
smiles= smiles[:3000]

# =============================================================================
# mols= [convert(x) for x in smiles]
# =============================================================================

m= metrics()
def f1(s):
    return m.calculateScore(convert(s))

def f2(s):
    return m.NP_score(convert(s))

def f3(s):
    return m.qed_score(convert(s))

def f4(s):
    if m.pains(convert(s)):
        return 1
    else: return 0


problem = Problem(num_of_variables=1,objectives=[f3, f1])
evolution = Evolution(problem)
evolution.evolve(fathers= smiles)



