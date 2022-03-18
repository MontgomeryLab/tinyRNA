# Feature Selection
We provide a Features Sheet (`features.csv`) in which you can define selection rules to more accurately capture counts for the small RNAs of interest. The parameters for these rules include attributes commonly used in the classification of small RNAs, such as length, strandedness, and 5' nucleotide. They are utilized at each alignment locus to determine which overlapping features should be assigned a portion of the read counts for the given sequence.

>**Important**: candidate features do not receive counts if they do not pass selection process described below

Selection occurs in three stages, with the output of each stage as input to the next:
1. Features are matched to rules based on their GFF3 column 9 attributes
2. Features which overlap an alignment are eliminated based on the hierarchy values and desired overlap characteristics defined in their mached rules
3. Remaining feature candidates are then selected based on the small RNA attributes of the alignment to which they are being assigned. These attributes are, again, defined in each feature's matched rules

## Stage 1: Feature Attribute Parameters
| features.csv columns: | Select for... | with value... |
|-----------------------|---------------| --- |

Each feature's column 9 attributes are searched for the key-value combinations defined in the `Select for...` and `with value...` columns. Features, and the rules they matched, are retained for overlap evaluation at alignment loci. Attribute keys are allowed to have multiple comma separated values, and these values are treated as a list; only one of the listed values needs to match the `with value...` to be considered a valid match to the rule.

For example, if a rule contained `Class` and `WAGO` in these columns, then a feature with attributes<br>`... ;Class=CSR,WAGO; ...` would be considered a match for the rule.

>**Tip**: These parameters are case sensitive. The capitalization in your rule must match the capitalization in your GFF3 files

## Stage 2: Hierarchy and Overlap Parameters
| features.csv columns: | Hierarchy | Match |
|-----------------------|-----------| --- |

This stage of selection is concerned with the interval overlap between alignments and features. **Overlap is determined in a strandless fashion.** See the `Strand` section in Stage 3 for refinement of selections by strand.

### Hierarchy
Each rule must be assigned a hierarchy value. This value is used to perform elimination when multiple features, or multiple feature-rule pairs, overlap an alignment locus.
- Each feature can have multiple hierarchy values if it matched more than one rule during Stage 1 selection
- Multiple rules are allowed to share the same value
- Only the lowest value is selected at each locus

>**Important:**
Let's take a step back. What exactly is the product of selection here? Not just a feature, but a feature _and_ a rule it had matched during Stage 1 selection. This is an important distinction because in Stage 3, only the **selected rule(s)** will be used to determine if the corresponding feature is an appropriate assignment based on the alignment's attributes.

You can use higher hierarchy values to exclude features that are not of interest.

>**Example:** suppose you have a miRNA locus embedded within a coding gene locus (within an intron for example). By assigning a hierarchy of 1 to miRNA and a hierarchy of 2 to coding genes, all small RNA counts from sequences matching to the miRNA would be excluded from total counts for the coding gene. Reversing the hierarchy such that miRNA had a hierarchy of 2 and coding genes had a hierarchy of 1 would instead exclude reads from sequences matching to the coding gene from total counts for the miRNA. If a hierarchy of 1 was assigned to both miRNAs and coding genes, counts for sequences matching both features would be split between them.

### Match
The match column allows you to specify which read alignments should be assigned based on how their start and end points overlap with candidate features. Candidates for each matched rule can be selected using the following options:
- `partial`: alignment overlaps feature by at least one base
- `full`: alignment does not extend beyond either terminus of the feature
- `exact`: alignment termini are equal to the feature's
- `5' anchored`: alignment's 5' end is equal to the corresponding terminus of the feature
- `3' anchored`: alignment's 3' end is equal to the corresponding terminus of the feature

The following diagrams demonstrate the strand semantics of these interval selectors. The first two options show separate illustrations for features on each strand for emphasis. All matches shown in the remaining three options apply to features on either strand.
![3'_anchored_5'_anchored](../images/3'_anchored_5'_anchored.png)
![Full_Exact_Partial](../images/Full_Exact_Partial.png)

## Stage 3: Alignment Attribute Parameters
| features.csv columns: | Strand | 5' End Nucleotide | Length |
|-----------------------|--------| --- | --- |

The final stage of selection is concerned with attributes of the alignment to which features are being assigned.

### Strand
- `sense`: the alignment strand must match the feature's strand for a match
- `antisense`: the alignment strand must not match the feature's strand for a match
- `both`: strand is not evaluated


### 5' End Nucleotide and Length
| Parameter | Single | List | Range | Wildcard |
| --- |:------:|:----:|:-----:|:--------:|
| 5' end nt |   ✓    |  ✓   |       |    ✓     |
| Length |    ✓    |   ✓   |   ✓    |     ✓     |

Examples:
- **Single**: `G` or `22`
- **List**: `C,G,U` or `25, 26` (spaces do not matter)
- **Range**: `20-25`
- **Wildcard**: `all`
- **Mixed**: `19, 21-23, 25-30`

>**Tip:** these parameters are **not** case sensitive.

>**Tip:** you may specify U and T bases in your rules. Uracil bases will be converted to thymine when your Features Sheet is loaded.

### Misc
| features.csv columns: | Alias by... | Feature Source |
| --- |-------------| --- |

You may specify an **Alias by...** which is a GFF3 column 9 attribute key you wish to represent each feature. The intention of this column is to provide a human-friendly name for each feature. The value associated with each feature's **Alias by...** attribute will be shown in the `Feature Name` column of the Feature Counts output table.  For example, if one of your rules specifies an alias of `sequence_name` and gene1's `sequence_name` attribute is "abc123", then gene1's `Feature Name` column in the Feature Counts table will read "abc123".

The **Feature Source** field of a rule is tied only to the **Alias by...**; rules are _not_ partitioned on a GFF file basis, and features parsed from these GFF files are similarly not partitioned as they all go into the same lookup table regardless of source. For each rule, aliases are built on a per-GFF file basis; that is, **Alias by...** values will only be gathered from their corresponding **Feature Source**. Additionally, each GFF file is parsed only once regardless of the number of times it occurs in the Features Sheet.

## Count Normalization
Small RNA reads passing selection will receive a normalized count increment. By default, read counts are normalized twice before being assigned to a feature. The second normalization step can be disabled in `run_config.yml` if desired. Counts for each small RNA sequence are divided: 
1. By the number of loci it aligns to in the genome.
2. By the number of _selected_ features for each of its alignments.

### The Details
You may encounter the following cases when you have more than one unique GFF file listed in your **Feature Source**s:
- If a feature is defined in one GFF3 file, then again in another but with differing attributes, rule and alias matches will be merged for the feature
- If a feature is defined in one GFF3 file, then again but under a different **Alias by...**, then both aliases are retained and treated as a list. All aliases will be present in the `Feature Name` column of the Feature Counts output table. They will be comma separated.

Discontinuous features and feature filtering support:
- Discontinuous features are supported (as defined by the `Parent` attribute key, or by a shared `ID` attribute value). Rule and alias matches of descendents are merged with the root parent's.
- Features can be filtered during GFF3 parsing by on their `source` and/or `type` columns, and these preferences can be specified in the Run Config file. These are inclusive filters. Only features matching the values specified will be retained for selection. An empty list allows all values.
- If a filtered feature breaks a feature lineage (that is, features chained via the `Parent` attribute), then the highest non-filtered ancestor is still considered to be the root parent. The lineage is maintained transparently but the filtered feature does not contribute to the domains of selection.