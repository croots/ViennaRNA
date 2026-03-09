#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <ViennaRNA/fold_compound.h>
#include <ViennaRNA/mfe.h>
#include <ViennaRNA/params/basic.h>
#include <ViennaRNA/utils/basic.h>
#include <ViennaRNA/eval/structures.h>
#include <ViennaRNA/subopt/wuchty.h>


static void
count_subopt_result(const char  *structure,
                    float       energy,
                    void        *data)
{
  if (structure)
    (*(int *)data)++;
}


#suite Subopt

#tcase DensityOfStates

/*
 * Regression test for the density_of_states out-of-bounds write.
 *
 * Root cause: in vrna_subopt_cb() with dangle_model 1 or 3, the MFE
 * structure is found with dangles temporarily forced to 2, then min_en is
 * re-evaluated under the original dangle model (1 or 3).  A suboptimal
 * structure enumerated under dangle=2 can score *lower* than that
 * re-evaluated min_en when rescored under dangle=1/3, giving
 *   e = (int)((structure_energy - min_en) * 10 - correction) < 0
 * Without the clamp "if (e < 0) e = 0;", density_of_states[e] increments
 * memory before the array, corrupting adjacent globals (e.g. a parameter-
 * file pointer) and crashing on the next ViennaRNA call.
 *
 * This sequence was found by fuzzing: with dangle=1 and delta=200
 * (2 kcal/mol), the structure "....................(((........)))." has
 * e = -2 before the clamp.
 */
#test test_subopt_dos_e_negative_clamped
{
  const char            *seq   = "GACUAGCUGUAAUCCGUAAAGCCGCCCAUGAGGCC";
  const int             delta  = 200;   /* 2.00 kcal/mol suboptimal window */
  vrna_md_t             md;
  vrna_fold_compound_t  *fc;
  int                   count1, count2;

  vrna_md_set_default(&md);
  md.dangles = 1;   /* triggers the dangle-model mismatch in subopt */

  /*
   * First call.  Without the clamp, this would write to
   * density_of_states[-2], corrupting adjacent memory.
   */
  fc      = vrna_fold_compound(seq, &md, VRNA_OPTION_MFE);
  count1  = 0;
  vrna_subopt_cb(fc, delta, count_subopt_result, &count1);
  vrna_fold_compound_free(fc);
  ck_assert_int_gt(count1, 0);

  /*
   * Second call on the same sequence.  If the first call corrupted memory,
   * this call would crash (e.g. when ViennaRNA tries to use a clobbered
   * pointer that was adjacent to density_of_states in BSS).  The result
   * must also match the first call.
   */
  fc      = vrna_fold_compound(seq, &md, VRNA_OPTION_MFE);
  count2  = 0;
  vrna_subopt_cb(fc, delta, count_subopt_result, &count2);
  vrna_fold_compound_free(fc);

  ck_assert_int_eq(count1, count2);
}


/*
 * Verify the pre-fix condition: at least one suboptimal structure for the
 * triggering sequence has structure_energy < min_en under dangle=1,
 * i.e. the unclamped index e would be negative.
 *
 * This test documents *why* the clamp is necessary without bypassing it.
 */
#test test_subopt_dos_e_negative_condition
{
  const char            *seq    = "GACUAGCUGUAAUCCGUAAAGCCGCCCAUGAGGCC";
  const int             delta   = 200;
  vrna_md_t             md;
  vrna_fold_compound_t  *fc;
  char                  *mfe_struc;
  double                min_en;
  float                 correction;
  int                   found_negative;
  int                   seq_len = (int)strlen(seq);

  vrna_md_set_default(&md);
  md.dangles = 1;
  fc = vrna_fold_compound(seq, &md, VRNA_OPTION_MFE);

  /*
   * Replicate the subopt.c logic:
   *   1. Force dangle=2, compute MFE structure.
   *   2. Restore dangle=1, re-evaluate to get min_en.
   */
  fc->params->model_details.dangles = 2;
  mfe_struc = (char *)vrna_alloc(sizeof(char) * (seq_len + 1));
  vrna_mfe(fc, mfe_struc);
  fc->params->model_details.dangles = 1;

  min_en      = vrna_eval_structure(fc, mfe_struc);
  correction  = (min_en < 0) ? -0.1f : 0.1f;
  free(mfe_struc);

  /*
   * Collect suboptimal structures (internally uses dangle=2 DP tables)
   * and rescore each under dangle=1.
   */
  vrna_subopt_solution_t *sols = vrna_subopt(fc, delta, 0, NULL);
  ck_assert_ptr_nonnull(sols);

  found_negative = 0;
  for (int i = 0; sols[i].structure != NULL; i++) {
    double  se    = vrna_eval_structure(fc, sols[i].structure);
    /* 10.0 converts kcal/mol to the integer deci-kcal units used as the
     * density_of_states index; mirrors the expression in vrna_subopt_cb(). */
    int     e_raw = (int)((se - min_en) * 10.0 - correction);
    if (e_raw < 0)
      found_negative = 1;
    free(sols[i].structure);
  }
  free(sols);
  vrna_fold_compound_free(fc);

  /* At least one structure must give e < 0 to confirm the clamp is needed */
  ck_assert_int_eq(found_negative, 1);
}


#main-pre
    srunner_set_tap(sr, "-");
