from pura.units import *
import pytest


@pytest.mark.parametrize(
    "unit_type_class,good_examples,failure_examples",
    [
        (Mass, ["25 g", "25 kg"], ["25 mol"]),
        (Volume, ["25 ml", "25 L"], ["25 g"]),
        (Temperature, [ureg.Quantity(25, "degC"), ureg.Quantity(25, "degF")], ["25 g"]),
        (Pressure, ["25 bar", "25 atm"], ["25 g"]),
        (Time, ["25 s", "25 min"], ["25 g"]),
        (Amount, ["25 mol", "25 mmol"], ["25 g"]),
        (MassFlow, ["25 g/s", "25 kg/min"], ["25 g"]),
        (VolumeFlow, ["25 ml/s", "25 L/min"], ["25 g"]),
        (MolarFlow, ["25 mol/s", "25 mmol/min"], ["25 g"]),
    ],
)
def test_pint_unit_model(unit_type_class, good_examples, failure_examples):
    class MockModel(PintModel):
        test_val: unit_type_class

    for pint_example in good_examples:
        m = MockModel(test_val=pint_example)
        assert m.test_val == ureg.Quantity(pint_example)
    for pint_example in failure_examples:
        with pytest.raises(ValueError):
            MockModel(test_val=pint_example)
