#include <silex/unit.hpp>

#include <cassert>
#include <iostream>
#include <vector>

namespace {
namespace sflint = silex::flint;

int run() {
    sflint::Fmpz radicand;
    fmpz_set_si(radicand.raw(), 2);

    silex::NumberField field =
            silex::NumberField::quadratic(sflint::FmpzConstRef(radicand));
    assert(field.is_defined());

    silex::Element theta(field);
    silex::Element epsilon(field);
    assert(theta.gen());
    assert(epsilon.add_si(theta, 1));

    silex::EmbeddingContext embeddings(field);
    assert(embeddings.is_defined());

    std::vector<silex::Element> units;
    units.emplace_back(field);
    assert(units[0].set(epsilon));

    sflint::ArbMat logs(1, 2);
    sflint::Arb regulator;
    const silex::ElementSpan unit_view(units.data(), units.size());
    assert(silex::unit_log_matrix(sflint::ArbMatRef(logs), embeddings,
                                  unit_view,
                                  silex::LogEmbeddingMode::product, 128));
    assert(silex::unit_regulator(sflint::ArbRef(regulator), embeddings,
                                 unit_view, 128));

    std::cout << "K = Q(sqrt(2))\n";
    std::cout << "unit = 1 + theta\n";
    std::cout << "rank = 1\n";
    std::cout << "log matrix rows = " << arb_mat_nrows(logs.raw())
              << ", cols = " << arb_mat_ncols(logs.raw()) << "\n";
    std::cout << "regulator enclosure: ";
    arb_printn(regulator.raw(), 20, 0);
    std::cout << "\n";
    return 0;
}

}  // namespace

int main() {
    return run();
}
