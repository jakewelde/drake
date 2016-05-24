This reeadme is to explain the figures in this directory.

Figures of the form drakeold_eg%02d.png were generated
using an older version of drake, one which sits under 
drc (eg. workstation commit f624742). Figures of the
form drakenew_eg%02.png were generated with a refactored
version of drake (master commit 6f47219).

One can clearly see huge improvements, especially to the
stability of the raycasting results. These images were
obtained by starting the simulated GelSight view frustum
far away, casting a ray at the object, then moving it
mostly or entirely that distance toward the object.
Thus, images where penetration occurs between the back
of the view frustum and the box (eg. drakeold_eg05.png)
are in error, and arise due to an instability in the 
raycasting algorithm.
