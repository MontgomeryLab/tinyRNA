# Feature Selection
Small RNAs can often be classified by sequence features, such as length, strandedness, and 5' nucleotide. We provide a Features Sheet (`features.csv`) in which you can define selection rules to more accurately capture counts for the small RNAs of interest.

Selection takes place for every feature associated with every alignment of every small RNA sequence. It occurs in two phases:
1. Against the candidate feature's attribute key/value pairs, as defined in your GFF3 file.
2. Against the alignment's small RNA attributes (strand relative to feature of interest, 5' end nucleotide, and length).

Each rule must be assigned a hierarchy value. A lower value indicates higher selection preference and multiple rules may share the same value. We utilize this value only during the first phase of selection; if multiple features match the attribute key/value pairs defined in your rules, then only the feature(s) with the lowest hierarchy values move to the second selection phase. The remaining features are discarded for the given read alignment. 

You can use higher hierarchy values to exclude features that are not of interest. For example, suppose you have a miRNA locus embedded within a coding gene locus (within an intron for example). By assigning a hierarchy of 1 to miRNA and a hierarchy of 2 to coding genes, all small RNA counts from sequences matching to the miRNA would be excluded from total counts for the coding gene. Reversing the hierarchy such that miRNA had a hierarchy of 2 and coding genes had a hierarchy of 1 would instead exclude reads from sequences matching to the coding gene from total counts for the miRNA. If a hierarchy of 1 was assigned to both miRNAs and coding genes, counts for sequences matching both features would be split between them.

>**Important**: candidate features do not receive counts if they do not pass selection process described below

## Feature Attribute Parameters
| Select for... | with value... | Hierarchy |
| --- | --- | --- |

The first round of selection is performed using this portion of the rule. Each rule must be assigned a hierarchy value, which applies only to this round of selection. A lower hierarchy value indicates higher selection preference and rules may share hierarchy values.

For each alignment, all associated features will first be examined to see if their column 9 attributes contain a "select for" key whose values contain a "with value". If multiple features match, then elimination will be performed using the hierarchy values of the rules they matched. The feature-rule pair(s) with the lowest hierarchy value will be selected for a second round of elimination. Internally, these matches are represented as `(hierarchy, rule, feature)` tuples which we call _hits_.

A feature may match multiple rules. When this happens, a _hit_ is produced for each matched rule and normal hierarchical elimination is performed. This means that the product of the first round of elimination is not just a list of features, but rather feature-rule pairs in the form of _hits_.

>**Tip**: These parameters are case sensitive. The capitalization in your rule must match the capitalization in your GFF files

## Read Attribute Parameters
| Strand | 5' End Nucleotide | Length | Match |
| --- | --- | --- | --- |

The second round of selection switches to attributes of the read to which features are being assigned. Recall that the first round of selection yields a list of _hits_. Now, only the rules contained in this list of _hits_ are used for selection. Contrast this with the first round in which the key-value pairs of all rules were considered.

### Strand
- `sense`: the alignment strand must match the feature's strand for a match
- `antisense`: the alignment strand must not match the feature's strand for a match
- `both`: strand is not evaluated

### Match
Valid values are `Partial` and `Full`. This parameter is referring to the required amount of overlap between a read alignment and a feature in order for that feature to be a candidate for selection. If a rule specifies a `Partial` Match value, then features which overlap a read alignment by at least one base will be considered for selection. If a rule specifies a `Full` Match value, then only features whose endpoints are fully contained by or equal to the read alignment will be considered for selection.

### 5' End Nucleotide and Length
| Parameter | Single | List | Range | Wildcard |
| --- | :---: | :---: | :---: | :---: |
| 5' end nt | X | X |  | X |
| Length | X | X | X | X |

Examples:
- **Single**: `G` or `22`
- **List**: `C,G,U` or `25, 26` (spaces do not matter)
- **Range**: `20-25`
- **Wildcard**: `all`
- **Mixed**: `19, 21-23, 25-30`

>**Tip:** these parameters are **not** case sensitive.

>**Tip:** you may specify U and T bases in your rules. Uracil bases will be converted to thymine when your Features Sheet is loaded.

### Misc
| Alias by... | Feature Source |
| --- | --- |

You may specify an **Alias by...** which is a GFF column 9 attribute key you wish to represent each feature. The intention of this column is to provide a human-friendly name for each feature. The value associated with each feature's **Alias by...** attribute will be shown in the `Feature Name` column of the Feature Counts output table.  For example, if one of your rules specifies an alias of `sequence_name` and gene1's `sequence_name` attribute is "abc123", then gene1's `Feature Name` column in the Feature Counts table will read "abc123".

The **Feature Source** field of a rule is tied only to the **Alias by...**; rules are _not_ partitioned on a GFF file basis, and features parsed from these GFF files are similarly not partitioned as they all go into the same lookup table regardless of source. For each rule, aliases are built on a per-GFF file basis; that is, **Alias by...** values will only be gathered from their corresponding **Feature Source**. Additionally, each GFF file is parsed only once regardless of the number of times it occurs in the Features Sheet.

## Count Normalization
Small RNA reads passing selection will receive a normalized count increment. By default, read counts are normalized twice before being assigned to a feature (these settings can be changed in `run_config.yml`). Counts for each small RNA sequence are divided: 
1. By the number of loci it aligns to in the genome.
2. By the number of selected features for each of its alignments.

### The Nitty Gritty Details
You may encounter the following cases when you have more than one unique GFF file listed in your **Feature Source**s:
- Attribute value lists are supported. Values which contain commas are treated as lists of values when parsing GFF files
- If a feature is defined in one GFF file, then again in another but with differing attributes, then those attribute values will be appended
- A feature may have multiple **Alias by...** attributes associated with it
- If a feature is defined in one GFF file, then again but under a different **Alias by...**, then both aliases are retained and treated as a list. All aliases will be present in the `Feature Name` column of the Feature Counts output table. They will be comma separated.