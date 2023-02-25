"""Contains classes for facilitating automated backward compatibility with older configuration files"""
import sys

from pkg_resources import resource_filename
from ruamel.yaml import CommentedMap, YAML
from importlib.metadata import version

TINYRNA_VERSION = version('tinyrna')


class YamlShop:
    """A collection of helper functions for editing YAML documents that
    have been loaded into memory with ruamel.yaml.

    At the time of writing, editing the document (CommentedMap) involves
    modifying three attributes in the document object:
      - A dictionary of the document's key/value pairs with original order.
        __get__, __set__, and other dictionary operations are passed to an
        OrderedDict base class which is extended to include insert()
      - ca: comment attribute data for each key
      - lc: line/column data for each key
    """

    @staticmethod
    def rename_key(doc: CommentedMap, key, new_name):
        if key not in doc: return
        YamlShop._rename_ordered_dict_key(doc, key, new_name)
        YamlShop._rename_linecol_key(doc, key, new_name)
        YamlShop._rename_comment_attribute_key(doc, key, new_name)

    @staticmethod
    def _rename_ordered_dict_key(doc: CommentedMap, key, new_name):
        target_idx = list(doc.keys()).index(key)
        value = doc[key]
        doc.pop(key, None)

        doc.insert(target_idx, new_name, value)

    @staticmethod
    def _rename_linecol_key(doc: CommentedMap, key, new_name):
        lc_value = doc.lc.data[key].copy()
        doc.lc.data.pop(key, None)

        doc.lc.data[new_name] = lc_value

    @staticmethod
    def _rename_comment_attribute_key(doc: CommentedMap, key, new_name):
        # Some keys don't have comment data
        if key not in doc.ca.items: return

        ca_value = doc.ca.items[key].copy()
        doc.ca.items.pop(key, None)
        doc.ca.items[new_name] = ca_value

    @staticmethod
    def add_mapping(doc: CommentedMap, prec_key, key_obj):
        key, val = tuple(*key_obj.items())
        key_order = list(doc.keys())

        if key in doc:
            return
        if prec_key not in doc:
            raise ValueError(f"Preceding_key not found: {prec_key}")

        target_idx = key_order.index(prec_key) + 1
        doc.insert(target_idx, key, val)
        key_order.insert(target_idx, key)

        if isinstance(val, (dict, list)):
            print("Adding complex mappings not yet supported.\n"
                  "Proceed with caution...", file=sys.stderr)

        # Comments & linebreaks are often (but not always!) attached to
        #   the preceding key. Move them down to the new key.
        inherit_prev = doc.ca.items[prec_key][2]
        doc.ca.items[key] = [None, None, inherit_prev, None]
        doc.ca.items[prec_key][2] = None


class RunConfigCompatibility:

    def __init__(self, config_obj: CommentedMap):
        self.config = config_obj.copy()
        self.vstart = (config_obj.get("version") or "0.0.0").strip("v")  # trust reported version for now
        definitions = resource_filename('tiny', 'templates') + "/compatibility/run_config_compatibility.yml"

        self.yaml = YAML()
        with open(definitions, 'r') as f:
            self.vmap = self.yaml.load(f)
            self.vorder = sorted(self.vmap.keys())

    def upgrade(self):
        itinerary = self.get_itinerary()
        if not itinerary: return self.config
        if any('block' in self.vmap[v] for v in itinerary):  # Todo: perform partial upgrades
            msg = "Your Run Config is out of date and couldn't be automatically updated.\n" \
                  "You will need to manually update it before you can proceed.\n" \
                  "Please see the release notes in the GitHub Releases section.\n" \
                  "\tRun Config version is: {}\n".format(self.vstart) + \
                  "\tInstalled tinyRNA version is: {}".format(TINYRNA_VERSION)
            raise ValueError(msg)

        for vers in itinerary:
            self.rename_keys(vers)
            # self.remove_keys(vers)    # Not yet implemented
            self.add_mappings(vers)

        self.config['version'] = latest = itinerary[-1]
        print(f"Automatically upgraded Run Config from {self.vstart} to {latest} in memory")
        return self.config

    def rename_keys(self, vers):
        for op in self.vmap[vers].get('rename', []):
            key, new_name = tuple(*op.items())
            YamlShop.rename_key(self.config, key, new_name)

    def remove_keys(self, vers):
        raise NotImplementedError()

    def add_mappings(self, vers):
        for op in self.vmap[vers].get('add', []):
            prec_key = op.pop('preceding_key')
            YamlShop.add_mapping(self.config, prec_key, op)

    def get_itinerary(self):
        try:
            # Document version is present in list
            starting_idx = self.vorder.index(self.vstart)
            return self.vorder[starting_idx + 1:]
        except ValueError:
            # Document version not present, infer position in timeline
            with_self = sorted([self.vstart, *self.vorder])
            starting_idx = with_self.index(self.vstart)
            return self.vorder[starting_idx:]

