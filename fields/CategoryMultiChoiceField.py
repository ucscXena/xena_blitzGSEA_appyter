from appyter.profiles.default.fields.MultiChoiceField import MultiChoiceField
from appyter.ext.json import try_json_loads

class CategoryMultiChoiceField(MultiChoiceField):
    @property
    def raw_value(self):
        value = try_json_loads(self.args['value'])
        if value == False:
            value = 'false'
        if not value:
            return []
        elif type(value) == list:
            return value
        elif type(value) == str:
            return [value]
        else:
            return [self.args['value']]

    def constraint(self):
        return True
        #return self.raw_value is not None and all(v in self.choices for v in self.raw_value)
