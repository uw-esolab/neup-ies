# Read before adding a JSON script!!!


The duct tape holding some of the code together can be a little brittle. It works! but requires strict syntax.

![json_structure][json_structure]

---
## **SSC_inputs**
---

As of Python 3.6, dictionary keys preserve order when retrieved. Therefore the ordering of entries in the `SSC_inputs` dictionary needs to be in a specific format.

We assume that there are data entries for each of three modules:
- a `Plant` module (`tcsmolten_salt`, `nuclear_tes`, etc.)
- a `Grid` module (typically `grid`)
- a `Financial` module (typically `single_owner`)

The **first** entry of `SSC_inputs` should be:
```
"compute_module_0" : <SSC_module_name>, 
```
for whatever module you are trying to simulate. All entries after that belong to the `Plant` module. 

The **next** entry (after the last `Plant` data entry is specified) should be:
```
"compute_module_1" : <SSC_grid_name>, 
```
followed by `Grid`-specific data entries.

Repeat the syntax for `Financial` module data entries.

**Finally**, the last entry of the `SSC_inputs` dictionary in the JSON script should be:
```
"number_compute_modules" : <N>, 
```
where *N* is the number of compute modules (typically 3).
This syntax ensures that we can run the `PySSCWrapper` correctly for debugging.

[json_structure]: json_structure.png