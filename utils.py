
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import sys, math, gzip, pickle
import os.path
import joblib
from rdkit.Chem.Descriptors import qed
from rdkit.Chem import FilterCatalog
from collections import defaultdict


class metrics():
    def __init__(self):
        self.fscore = joblib.load('np/publicnp.model')
        data = joblib.load('sa/myZip.pkl')
        outDict = {}
        for i in data:
            for j in range(1, len(i)):
                outDict[i[j]] = float(i[0])
                
        self._fscores = outDict
        
    def qed_score(self, mol):
        return qed(mol)
    
    def pains(self, mol):
        """
        Baell and Holloway (2010) New Substructure Filters for Removal of Pan
        Assay Interference Compounds (PAINS) from Screening Libraries and for
        Their Exclusion in Bioassays
        This filter finds promiscuous compounds that are likely to show activity
        regardless of the target.
        Returns:
            Boolean of whether the molecule triggers the PAINS filter.
        """
        params = FilterCatalog.FilterCatalogParams()
        params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
        catalog = FilterCatalog.FilterCatalog(params)
        return catalog.HasMatch(mol)
    
    def NP_score(self, mol):

        if mol is None:
                raise ValueError('invalid molecule')
        fp = Chem.rdMolDescriptors.GetMorganFingerprint(mol, 2)
        bits = fp.GetNonzeroElements()
    
        # calculating the score
        score = 0.
        for bit in bits:
            score += self.fscore.get(bit, 0)
        score /= float(mol.GetNumAtoms())
    
        # preventing score explosion for exotic molecules
        if score > 4:
            score = 4. + math.log10(score - 4. + 1.)
        if score < -4:
            score = -4. - math.log10(-4. - score + 1.)
        return score     
        
    def numBridgeheadsAndSpiro(self, mol, ri=None):
        nSpiro = rdMolDescriptors.CalcNumSpiroAtoms(mol)
        nBridgehead = rdMolDescriptors.CalcNumBridgeheadAtoms(mol)
        return nBridgehead, nSpiro
    
    
    def calculateScore(self, m):
    
        # fragment score
        fp = rdMolDescriptors.GetMorganFingerprint(m, 2)  # <- 2 is the *radius* of the circular fingerprint
        fps = fp.GetNonzeroElements()
        score1 = 0.
        nf = 0
        for bitId, v in fps.items():
            nf += v
            sfp = bitId
            score1 += self._fscores.get(sfp, -4) * v
        score1 /= nf
    
        # features score
        nAtoms = m.GetNumAtoms()
        nChiralCenters = len(Chem.FindMolChiralCenters(m, includeUnassigned=True))
        ri = m.GetRingInfo()
        nBridgeheads, nSpiro = self.numBridgeheadsAndSpiro(m, ri)
        nMacrocycles = 0
        for x in ri.AtomRings():
            if len(x) > 8:
                nMacrocycles += 1
    
        sizePenalty = nAtoms**1.005 - nAtoms
        stereoPenalty = math.log10(nChiralCenters + 1)
        spiroPenalty = math.log10(nSpiro + 1)
        bridgePenalty = math.log10(nBridgeheads + 1)
        macrocyclePenalty = 0.
        # ---------------------------------------
        # This differs from the paper, which defines:
        #  macrocyclePenalty = math.log10(nMacrocycles+1)
        # This form generates better results when 2 or more macrocycles are present
        if nMacrocycles > 0:
            macrocyclePenalty = math.log10(2)
    
        score2 = 0. - sizePenalty - stereoPenalty - spiroPenalty - bridgePenalty - macrocyclePenalty
    
        # correction for the fingerprint density
        # not in the original publication, added in version 1.1
        # to make highly symmetrical molecules easier to synthetise
        score3 = 0.
        if nAtoms > len(fps):
            score3 = math.log(float(nAtoms) / len(fps)) * .5
    
        sascore = score1 + score2 + score3
    
        # need to transform "raw" value into scale between 1 and 10
        min = -4.0
        max = 2.5
        sascore = 11. - (sascore - min + 1) / (max - min) * 9.
        # smooth the 10-end
        if sascore > 8.:
            sascore = 8. + math.log(sascore + 1. - 9.)
        if sascore > 10.:
            sascore = 10.0
        elif sascore < 1.:
            sascore = 1.0
    
        return sascore


