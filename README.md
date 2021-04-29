# Search


Provides the types for objects and scenario management. This module also includes plot recipes for the various object types.

## Types
- `SimObject`: the abstract type that defines a basic interface for all objects within the sim environment.
    - `Target`: the subtype of `SimObject` that specializes for the search targets.
    - `Agent`: the subtype of `SimObject` that specializes for the searching agent.
    - `Sensor`: the subtype of `SimObject` that specializes for search sensors.
- `Scenario`: the type that contains an `Agent`, `Target`, `Sensors`, and other parameters.