use bitmasks
BEGIN_SHELL [ /usr/bin/env python ]
from generate_h_apply import *

s = H_apply("just_1h_1p")
s.set_selection_pt2("epstein_nesbet_2x2")
s.unset_skip()
s.filter_only_1h1p()
print s


s = H_apply("all_but_1h_and_1p")
s.set_selection_pt2("epstein_nesbet_2x2")
s.unset_skip()
s.filter_1h()
s.filter_1p()
print s



s = H_apply("standard")
s.set_selection_pt2("epstein_nesbet")
s.unset_skip()
print s

s = H_apply("just_mono")
s.set_selection_pt2("epstein_nesbet_2x2")
s.unset_skip()
s.unset_double_excitations()
print s



s = H_apply("just_mono_no_1h_no_1p")
s.set_selection_pt2("epstein_nesbet_2x2")
s.unset_skip()
s.unset_double_excitations()
s.filter_1h()
s.filter_1p()
print s

s = H_apply("just_mono_no_1h_no_1p_no_2p")
s.set_selection_pt2("epstein_nesbet_2x2")
s.unset_skip()
s.unset_double_excitations()
s.filter_1h()
s.filter_1p()
s.filter_2p()
print s


s = H_apply("h_core_just_mono")
s.set_selection_pt2("h_core")
s.unset_skip()
s.unset_double_excitations()
print s


s = H_apply("just_mono_no_double_hole_pt2")
s.set_perturbation("epstein_nesbet_2x2")
s.unset_double_excitations()
s.set_filter_holes()
print s


END_SHELL

