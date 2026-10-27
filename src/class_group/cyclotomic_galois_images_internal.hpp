#pragma once

#include <silex/factor_base.hpp>
#include <silex/flint/fmpz_mat.hpp>
#include <silex/number_field.hpp>
#include <silex/order.hpp>

#include <vector>

namespace silex::detail {

// Galois images of class-group relations on cyclotomic fields.
//
// Sources.  Reference buch2.c (2.17.3), add_rel (lines 2380-2410) tries the
// image of every relation for which add_rel_i (lines 2284-2376) returns
// k > 0 under each non-trivial automorphism, and passes each image to
// add_rel_i with m = NULL.  add_rel_i returns k > 0 for a relation that
// raises the modular rank, for a dependent relation stored under the
// surplus allowance relsup (lines 2351-2361), and for every new relation
// once the rank is full (cache->missing == 0, line 2350); so the reference
// keeps every image that is not already known (already_known, line 2294)
// once the rank is full.  Reference Rel_add.jl:class_group_add_relation
// (v0.38.6, lines 97-131) adds images only for a relation that add_gen!
// reports as a new generator, and keeps an image only when it is again
// one, so it adds none after full rank.
//
// Silex adds images without a cap on fields of degree at least
// kLargeFieldRelationSearchMinDegree (30, class_group_storage_internal.hpp),
// the crossover of the one-random-ideal-per-pass rule.  Below it, only a
// field defined by exactly Phi_n (not psi_n) adds images, only when
// capped_galois_images_selected admits it (degree at most 20, or a
// factor-base bound of at least 300), and only under a capped rule: within
// one relation round of the LLL route, kernel images (images stored at full
// modular rank) stop once the round's relation target plus |G| - 1 of them
// are stored, and those kernel images count toward the kernel-unit target.
// The crossover, the gate and the cap are Silex's own rules, set by
// measurement (the references add images at every degree, and the
// reference relation cache, buch2.c:2295 with pre_allocate, lines 351-355,
// is never this tight): uncapped images stored after full rank cost more
// relation work than they save on several cyclotomic corpus rows below
// degree 30.
//
// Silex follows buch2.c here.  It tries the images of every stored
// relation for which add_rel_i would return k > 0, and keeps an image by
// the rule of add_rel_i with m = NULL: before full modular rank only when
// it raises the rank (the generator path has no surplus allowance and
// skips dependent rows; the admission cache, the port of the reference
// relation cache, also applies relsup), and at full rank whenever it is
// not already known.  The bookkeeping on the generator path is Silex's
// own, not the reference's: the images of one relation are stored as a
// batch: they run no presentation publication of their
// own, they are left out of the relation cap max_relations, and, from
// degree 30 on, images stored without raising the rank are left out of the
// kernel-unit target (ClassGroupContextStorage::
// galois_image_kernel_relations; under the capped rule below degree 30
// they count toward it and are counted in
// galois_image_capped_kernel_relations instead).  The
// exclusion holds only where the caps and targets use
// ClassGroupRelationSearchAccess::capped_relation_count and the kernel
// counts that subtract the image counters (the LLL route and the goal
// checks of relation_search.cpp); other readers of relation_count() and
// relation_kernel_unit_count() (the relation-cap stops of the search
// driver, the completion scheduler, relation_kernel_target_open) still
// include images.  Images stored on the admission-cache path are counted
// in the same counters, so they too are left out of that later cap and
// target.  On a grh
// class/unit request for a field defined by Phi_n the cyclotomic units are
// seeded into the unit context before relation search (source
// add_cyclotomic_units, buch2.c:3626-3636; src/order_unit/compute.cpp)
// and count toward the kernel-unit need there, as the source counts them
// in zc (buch2.c:4060).  On a field defined by psi_n the real cyclotomic
// units (src/unit/cyclotomic_units.cpp) are seeded there in the same way;
// the source seeds nothing on such fields.
//
// FB_aut_perm (buch2.c lines 197-249) maps every factor-base prime to a
// factor-base prime above the same rational prime with the same residue
// degree, so the reference requires a factor base stable under the
// automorphisms.  Both references take the automorphisms from a general
// automorphism routine (buch2.c:automorphism_matrices, line 2724, through
// galoisconj).  Silex has no general automorphism algorithm, so this route
// covers only a field whose defining polynomial is exactly the cyclotomic
// polynomial Phi_n, or exactly psi_n, the minimal polynomial of
// 2 cos(2 pi / n).  On Phi_n, K = Q(zeta_n) with zeta_n = x, and Gal(K/Q)
// is (Z/nZ)^* acting by sigma_k(x) = x^k for k coprime to n (Washington,
// Introduction to Cyclotomic Fields, 2nd ed., Theorem 2.5).  On psi_n,
// K = Q(zeta_n)^+ with x = zeta_n + zeta_n^-1, Gal(K/Q) is the quotient
// (Z/nZ)^* / {+-1}, and the restriction of sigma_k maps x to
// zeta_n^k + zeta_n^-k = D_k(x), the Dickson polynomial
// (D_0 = 2, D_1 = x, D_(j+1) = x D_j - D_(j-1)), for 1 < k < n / 2.
//
// Silex does not permute exponent vectors as the references do.  It forms
// the image element sigma_k(alpha) exactly (its order coordinates times the
// exact integral image matrix of sigma_k) and factors it through the
// ordinary exact relation path, so every image relation is checked exactly
// like any searched relation, and a relation whose image is not smooth over
// the base is simply not added.  The reference stops the image loop of a
// relation at the first image row that is already known (Rel_add.jl:122);
// Silex tries every image, so it may add independent images the reference
// would skip.  Each is still an exact relation.
//
// Factor-base stability.  A field automorphism maps a prime P above p to a
// prime above the same p with the same residue degree, and it permutes the
// primes above p (Neukirch, Algebraic Number Theory, Proposition I.9.1);
// moreover v_{sigma(P)}(sigma(alpha)) = v_P(alpha).  So if the factor base
// contains every prime ideal above each rational prime it uses, it is stable
// under every automorphism, and the image of a relation that is smooth over
// the base is again smooth over it.  Silex records exactly that property per
// rational-prime block (FactorBaseBlockAccess::rational_prime_block_is_complete);
// the route runs only when every block is complete.
struct CyclotomicGaloisImages {
    // `order` is the order whose coordinates the matrices act on; the data
    // is rebuilt when a different order is passed.
    bool initialized = false;
    bool available = false;
    // True when the field is defined by psi_n (the maximal real subfield of
    // Q(zeta_n)), false for Phi_n.
    bool real_subfield = false;
    ulong conductor = 0;
    Order order;
    // image_matrices[i] is the integral matrix of sigma_{exponents[i]} on the
    // basis of `order`: row j holds the coordinates of sigma(b_j), so the
    // coordinate row c of alpha maps to c * image_matrices[i].
    std::vector<ulong> exponents;
    std::vector<flint::FmpzMat> image_matrices;
    // Counters of image candidates tried and image relations appended.
    slong image_candidates = 0;
    slong image_relations = 0;

    void clear() noexcept {
        initialized = false;
        available = false;
        real_subfield = false;
        conductor = 0;
        order.clear();
        exponents.clear();
        image_matrices.clear();
        image_candidates = 0;
        image_relations = 0;
    }
};

// Sets `n` to the index of the cyclotomic polynomial equal to the defining
// polynomial of `field`, or to 0 when the defining polynomial is not exactly
// some Phi_n.  The comparison is exact (FLINT fmpz_poly_is_cyclotomic, then
// equality with fmpz_poly_cyclotomic(n)).  Fails only on an undefined field.
bool cyclotomic_field_index(ulong& n, const NumberField& field) noexcept;

// Builds the automorphisms sigma_k : x -> x^k, k in [2, n) coprime to n, of a
// field of degree at least 4 defined by exactly Phi_n, or
// sigma_k : x -> D_k(x), 2 <= k < n / 2 coprime to n, of one defined by
// exactly psi_n (real_subfield), as exact FieldHom generator images
// (FieldHom::set_generator_image checks f(image) = 0 in K for the defining
// polynomial f), and their integral image matrices on the basis of `order`
// (OrderHom::image_matrix, which checks integrality).  Leaves `available`
// false for any other field.  Degree-two cyclotomic fields are imaginary
// quadratic: their class-group route keeps one prime per split rational
// prime (an unstable base), so they are excluded.  Fails only when an exact
// construction step fails on a field that is exactly Phi_n or psi_n.
bool prepare_cyclotomic_galois_images(CyclotomicGaloisImages& out,
                                      const Order& order) noexcept;

// True when every rational-prime block of `base` contains all prime ideals
// above its rational prime; see the stability note above.
bool factor_base_rational_prime_blocks_complete(
        bool& complete,
        const FactorBase& base) noexcept;

// out = coordinates * images.image_matrices[index].
bool cyclotomic_galois_image_coordinates(
        flint::FmpzMat& out,
        const CyclotomicGaloisImages& images,
        slong index,
        flint::FmpzMatConstRef coordinates) noexcept;

}  // namespace silex::detail
