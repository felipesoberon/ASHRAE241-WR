# ASHRAE 241 ECAi and ASHRAE 62.1 / 170 outdoor air: preliminary comparison

This working comparison places the published ASHRAE 241 ECAi values used in
the model beside candidate ASHRAE 62.1-2025 ventilation-rate procedure (VRP)
breathing-zone outdoor-air values and representative ASHRAE 170-2025 outdoor-air
rates and the baseline outdoor-air rates published by Jones et al. All are
expressed in L/s per person, but outdoor air and ECAi serve different purposes
and must not be treated as equivalent requirements.

## Method

For each selected 62.1 occupancy, calculate `Vbz / Pz = Rp + Ra × 100 / D`,
where `Rp` is in L/s per person, `Ra` in L/s per m², and `D` is the table's
default number of occupants per 100 m². Values below use the SI entries in
Table 6-1 and are rounded to two decimal places. No particular floor area is
needed: it cancels when the table's default density is used.

These are breathing-zone outdoor-air rates. Actual zone or system outdoor-air
requirements can differ with design population, area, air-distribution
effectiveness, system configuration and applicable provisions. They are not
total equivalent clean air including recirculation filtration or other cleaning.

For the representative healthcare cases, convert 170 minimum outdoor air changes
per hour (ACH) using the 241 model's volume `V` (m³) and occupancy `P`:
`outdoor air / person = outdoor ACH × V × 1000 / (3600 × P)`. Standard 170
does not itself provide these model occupancies or volumes. The article's
mapping is many-to-one: the chosen 170 row is one example, not a unique rate
for every real space in that 241 category. Total supply-air requirements,
filtration and other provisions must be assessed separately.

## Working comparison

| ASHRAE 241 model space | 241 ECAi, L/s/person | Candidate category and outdoor-air source | Our outdoor air, L/s/person | Jones baseline, L/s/person |
|---|---:|---|---:|---:|
| Cell | 15 | Cell | 4.90 | 4.9 |
| Dayroom | 20 | Dayroom | 3.50 | 3.5 |
| Food and beverage facilities | 30 | Restaurant dining rooms (candidate) | 5.14 | 5.6 |
| Gym | 40 | Gym, sports arena play area (candidate) | 22.86 | 11 |
| Office | 15 | Office space | 8.50 | 8.5 |
| Retail | 20 | Sales, except listed exceptions (candidate) | 7.80 | 7.8 |
| Transportation waiting | 30 | Transportation waiting | 4.25 | 6.8 |
| Classroom | 20 | Classroom, age 9+ (candidate) | 6.71 | 7.4 |
| Lecture hall | 25 | Lecture hall, fixed seats | 4.01 | 2.9 |
| Manufacturing | 25 | Manufacturing, hazardous materials not used (candidate) | 17.86 | 9.3 |
| Sorting, packing, light assembly | 10 | Sorting, packing, light assembly | 12.37 | 11 |
| Warehouse | 10 | 62.1 warehouse, using assumed 6 m ceiling | 8.00 | 6.5 |
| Exam room | 20 | 170 inpatient general examination room, Table 7-1 | 7.59 | 7.5 |
| Group treatment area | 35 | 170 outpatient cancer treatment area, Table 8-1 | 7.50 | 7.5 |
| Patient room | 35 | 170 inpatient patient room, Table 7-1 | 15.00 | 15 |
| Resident room | 25 | 170 nursing-home resident room, Table 9-1 | 15.00 | 15 |
| Waiting room | 45 | 170 emergency-department public waiting area, Table 8-1 | 5.00 | 5 |
| Auditorium | 25 | Auditorium seating area | 2.71 | 3.5 |
| Place of religious worship | 25 | Places of religious worship | 2.75 | 2.8 |
| Museum | 30 | Museums/galleries | 4.71 | 4.6 |
| Convention | 30 | 62.1 gambling casinos (representative) | 4.79 | 4.6 |
| Spectator area | 25 | Spectator areas | 4.01 | 4.7 |
| Lobbies | 25 | Public assembly lobbies (candidate) | 2.70 | 3.5 |
| Common space | 25 | 62.1 hotel lobbies/prefunction (representative) | 4.80 | 1 |
| Dwelling unit | 15 | 62.1 hotel bedroom/living room (representative) | 5.50 | 6.8 |

The non-healthcare values are 62.1-2025 breathing-zone rates at its default
occupant density, except for the warehouse assumption below. The healthcare values are 170-2025 outdoor ACH converted
using 241 model inputs. Their shared units do not make the methods equivalent.

The Jones column reproduces the baseline outdoor-air rates in Jones et al.,
Table 4, as published; Table 2 provides the underlying model inputs. The
contour plots use these published baseline values as their outdoor-air
reference, including Standard 62.2 for residential cases, alongside the
Standard 241 ECAi reference. The earlier calculations remain here as a
separate comparison based on the stated 2025 mappings and assumptions.

For common space, Jones Table 4 reports 1.0 L/s per person, whereas its
Table 2 inputs (2.5 L/s per person, 0.3 L/s per m², 500 m² and 150
occupants) imply 3.5 L/s per person. The plotted reference uses the
published Table 4 value pending clarification.

## Warehouse assumption

Standard 62.1-2025 Table 6-1 specifies a people outdoor-air rate of 5 L/s per
person and an area outdoor-air rate of 0.3 L/s per m² for warehouses, but does
not specify a default occupant density. For this example, use the 241 model's
volume of 1,200 m³ and 20 occupants, and assume a 6 m ceiling height. The
resulting estimated floor area is 1,200 / 6 = 200 m². The breathing-zone
outdoor airflow is 5 × 20 + 0.3 × 200 = 160 L/s, or **8.00 L/s per person**.
The ceiling height and resulting area are assumptions, not values specified by
the standards; an actual warehouse's geometry and occupancy can change the rate.

## Healthcare assumptions

| Healthcare model space | 170 outdoor ACH | Model volume, m³ | Model occupants | Computed L/s/person |
|---|---:|---:|---:|---:|
| Exam room | 2 | 41 | 3 | 7.59 |
| Group treatment area | 2 | 270 | 20 | 7.50 |
| Patient room | 2 | 81 | 3 | 15.00 |
| Resident room | 2 | 81 | 3 | 15.00 |
| Waiting room | 2 | 270 | 30 | 5.00 |

For example, the selected emergency-department waiting area also has a total
supply-air minimum, which is not represented by the outdoor-air-only value.

Candidate matches need review against the actual space use and geometry. For
example, 62.1 also lists classrooms for ages 5–8 at **7.40 L/s per person**,
and an office main-entry lobby has a different default density from the public
assembly lobby selected above. The category name alone does not establish the
correct match. For healthcare spaces within its scope, ASHRAE 170 rates apply;
62.1 Appendix E addresses certain outpatient facilities outside that scope.

The 241 model's own occupant counts and geometry need not match the 62.1
default density. Any contour for a specific modeled space should therefore
use a documented crosswalk, actual floor area and population, and applicable
system assumptions rather than automatically substituting these defaults.
The convention, common-space and dwelling-unit rows are illustrative article
mappings; they do not establish the correct prerequisite standard or reference
rate for every convention, residential common space or dwelling unit. In
particular, separately check whether Standard 62.2 applies to a dwelling.

## Sources

- ASHRAE 62.1-2025, Section 6.2.1.1 and Table 6-1; Appendix E for applicable outpatient facilities.
- ASHRAE/ASHE 170-2025, Tables 7-1, 8-1 and 9-1.
- McNulty, English and Zaatari, “Mapping Space Types Across ASHRAE Standards 62.1, 170, & 241,” ASHRAE Journal, June 2025, Table 2 and supplementary 62.1 mapping. The mapping uses 62.1-2022; rates above are independently read from 62.1-2025.
- ASHRAE 241 published comparator values as recorded in [`../../src/category_data.cpp`](../../src/category_data.cpp). These are the model's reference values, not values derived by this comparison.
- Jones et al., “Risk modeling for ASHRAE Standard 241-2023 — Control of infectious aerosols,” Building and Environment 283 (2025), 113318, Tables 2 and 4. https://doi.org/10.1016/j.buildenv.2025.113318
