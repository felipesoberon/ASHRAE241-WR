# ASHRAE 241 ECAi and ASHRAE 62.1-2025 outdoor air: preliminary comparison

This working comparison places the published ASHRAE 241 ECAi values used in
the model beside candidate ASHRAE 62.1-2025 ventilation-rate procedure (VRP)
breathing-zone outdoor-air values. Both are expressed in L/s per person, but
they serve different purposes and must not be treated as equivalent requirements.

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

## Working comparison

| ASHRAE 241 model space | 241 ECAi, L/s/person | Candidate 62.1-2025 category | 62.1 Vbz/Pz, L/s/person |
|---|---:|---|---:|
| Cell | 15 | Cell | 4.90 |
| Dayroom | 20 | Dayroom | 3.50 |
| Food and beverage facilities | 30 | Restaurant dining rooms (candidate) | 5.14 |
| Gym | 40 | Gym, sports arena play area (candidate) | 22.86 |
| Office | 15 | Office space | 8.50 |
| Retail | 20 | Sales, except listed exceptions (candidate) | 7.80 |
| Transportation waiting | 30 | Transportation waiting | 4.25 |
| Classroom | 20 | Classroom, age 9+ (candidate) | 6.71 |
| Lecture hall | 25 | Lecture hall, fixed seats | 4.01 |
| Manufacturing | 25 | Manufacturing, hazardous materials not used (candidate) | 17.86 |
| Sorting, packing, light assembly | 10 | Sorting, packing, light assembly | 12.37 |
| Warehouse | 10 | Warehouses | Not calculated: no default density given |
| Exam room | 20 | Healthcare: check scope and ASHRAE 170 | Not assigned |
| Group treatment area | 35 | Healthcare: check scope and ASHRAE 170 | Not assigned |
| Patient room | 35 | Healthcare: check scope and ASHRAE 170 | Not assigned |
| Resident room | 25 | Healthcare: check scope and ASHRAE 170 | Not assigned |
| Waiting room | 45 | Healthcare: check scope and ASHRAE 170 | Not assigned |
| Auditorium | 25 | Auditorium seating area | 2.71 |
| Place of religious worship | 25 | Places of religious worship | 2.75 |
| Museum | 30 | Museums/galleries | 4.71 |
| Convention | 30 | No unique match selected | Not assigned |
| Spectator area | 25 | Spectator areas | 4.01 |
| Lobbies | 25 | Public assembly lobbies (candidate) | 2.70 |
| Common space | 25 | No unique match selected | Not assigned |
| Dwelling unit | 15 | Check residential scope and ASHRAE 62.2 | Not assigned |

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

## Sources

- ASHRAE 62.1-2025, Section 6.2.1.1 and Table 6-1; Appendix E for applicable outpatient facilities. Local document: `Standards Library - Documents/ASHRAE/62.1/01 Standards/2025 - ASHRAE - 62.1 Ventilation and Acceptable Indoor Air Quality.pdf`.
- ASHRAE 241 published comparator values as recorded in [`../../src/category_data.cpp`](../../src/category_data.cpp). These are the model's reference values, not values derived by this comparison.
