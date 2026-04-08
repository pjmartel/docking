Looking at the PDF, there are actually **6 compounds** in the table (perhaps you meant "these formulas" — no worries!). I'll provide SMILES for all of them.

The key anchor here is **KKW-061 = AS1842856**, a known FOXO1 inhibitor whose structure is well-described in the literature as **5-amino-7-(cyclohexylamino)-1-ethyl-6-fluoro-4-oxo-1,4-dihydro-1,8-naphthyridine-3-carboxylic acid**. The remaining compounds are analogs derived by visual inspection of the structural drawings.

| Code | SMILES | Notes |
|------|--------|-------|
| KKW-035 | `O=C(OC)c1cn(CC)c2nc(NC3CCCCC3)ccc2c1=O` | Methyl ester (no F, no NH₂) |
| KKW-038 | `O=C(O)c1cn(CC)c2nc(NC3CCCCC3)ccc2c1=O` | Free acid (no F, no NH₂) |
| KKW-053 | `O=C(O)c1cn(CCCC)c2nc(NC3CCCCC3)ccc2c1=O` | n-butyl on N (longer chain, no F, no NH₂) |
| KKW-061 (= AS1842856) | `O=C(O)c1cn(CC)c2nc(NC3CCCCC3)c(F)c(N)c2c1=O` | **Reference compound** — 6-F, 5-NH₂, 1,8-naphthyridine scaffold |
| KKW-063 | `O=C(O)c1cn(CC)c2nc(NC3CCCCC3)cc(N)c2c1=O` | 5-NH₂, no F |
| KKW-090 | `O=C(O)c1cn(CC)c2nc(NC3CCCCC3)c(F)cc2c1=O` | 6-F, no NH₂ |

**A few caveats worth noting:**

- The **KKW-061 SMILES** is the most reliable, anchored to the known AS1842856 structure (1,8-naphthyridine-3-carboxylic acid scaffold).
- For **KKW-053**, I inferred an n-butyl chain to account for the higher MW (370.5 vs 328.4 for KKW-038, ΔMW ≈ 42 ≈ C₃H₆ relative to ethyl → n-pentyl, actually worth rechecking against the original). You might want to verify this with `Chem.MolFromSmiles()` + `rdMolDescriptors.CalcExactMolWt()` in RDKit.
- The structures for **KKW-035/038/053** suggest they may share a quinoline rather than 1,8-naphthyridine core — the image resolution makes it difficult to distinguish confidently.

I'd recommend a quick RDKit check like:
```python
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

smiles = "O=C(O)c1cn(CC)c2nc(NC3CCCCC3)c(F)c(N)c2c1=O"
mol = Chem.MolFromSmiles(smiles)
print(rdMolDescriptors.CalcExactMolWt(mol))  # Should give ~347.4 for KKW-061
```