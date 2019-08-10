// Compute the exp(1.1 * x) / 12 for x real.
const e_float one_pt_one(1.1);
const e_float exp_x_times_one_pt_one = ef::exp(x * one_pt_one);

const e_float y = exp_x_times_one_pt_one / static_cast<INT32>(12);

