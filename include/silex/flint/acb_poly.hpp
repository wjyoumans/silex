#pragma once

#include <flint/acb_poly.h>

#include <utility>

namespace silex::flint {

class AcbPolyRef;
class AcbPolyConstRef;

class AcbPoly {
public:
    AcbPoly() noexcept { acb_poly_init(value_); }
    ~AcbPoly() noexcept { acb_poly_clear(value_); }

    AcbPoly(const AcbPoly&) = delete;
    AcbPoly& operator=(const AcbPoly&) = delete;

    AcbPoly(AcbPoly&& other) noexcept {
        acb_poly_init(value_);
        acb_poly_swap(value_, other.value_);
    }

    AcbPoly& operator=(AcbPoly&& other) noexcept {
        if (this != &other) {
            acb_poly_swap(value_, other.value_);
        }
        return *this;
    }

    void swap(AcbPoly& other) noexcept { acb_poly_swap(value_, other.value_); }
    acb_poly_t& raw() noexcept { return value_; }
    const acb_poly_t& raw() const noexcept { return value_; }

private:
    acb_poly_t value_;
};

class AcbPolyConstRef {
public:
    explicit AcbPolyConstRef(const acb_poly_t& value) noexcept : value_(value) {}
    explicit AcbPolyConstRef(const AcbPoly& value) noexcept : value_(value.raw()) {}
    const acb_poly_struct* raw() const noexcept { return value_; }

private:
    const acb_poly_struct* value_;
};

class AcbPolyRef {
public:
    explicit AcbPolyRef(acb_poly_t& value) noexcept : value_(value) {}
    explicit AcbPolyRef(AcbPoly& value) noexcept : value_(value.raw()) {}
    acb_poly_struct* raw() noexcept { return value_; }
    const acb_poly_struct* raw() const noexcept { return value_; }

private:
    acb_poly_struct* value_;
};

class AcbPolyEvaluationTree {
public:
    AcbPolyEvaluationTree() noexcept = default;

    ~AcbPolyEvaluationTree() noexcept { clear(); }

    AcbPolyEvaluationTree(const AcbPolyEvaluationTree&) = delete;
    AcbPolyEvaluationTree& operator=(const AcbPolyEvaluationTree&) = delete;

    AcbPolyEvaluationTree(AcbPolyEvaluationTree&& other) noexcept {
        swap(other);
    }

    AcbPolyEvaluationTree& operator=(
            AcbPolyEvaluationTree&& other) noexcept {
        if (this != &other) {
            AcbPolyEvaluationTree tmp(std::move(other));
            swap(tmp);
        }
        return *this;
    }

    void swap(AcbPolyEvaluationTree& other) noexcept {
        std::swap(tree_, other.tree_);
        std::swap(length_, other.length_);
        std::swap(precision_, other.precision_);
    }

    void clear() noexcept {
        if (tree_ != nullptr) {
            _acb_poly_tree_free(tree_, length_);
        }
        tree_ = nullptr;
        length_ = 0;
        precision_ = 0;
    }

    bool build(acb_srcptr points,
               slong length,
               slong precision) noexcept {
        if (points == nullptr || length <= 0 || precision <= 0) {
            return false;
        }

        AcbPolyEvaluationTree candidate;
        candidate.tree_ = _acb_poly_tree_alloc(length);
        if (candidate.tree_ == nullptr) {
            return false;
        }
        _acb_poly_tree_build(candidate.tree_, points, length, precision);
        candidate.length_ = length;
        candidate.precision_ = precision;
        swap(candidate);
        return true;
    }

    bool evaluate(acb_ptr out,
                  AcbPolyConstRef polynomial,
                  slong precision) const noexcept {
        if (out == nullptr || tree_ == nullptr || length_ <= 0 ||
            precision != precision_) {
            return false;
        }
        _acb_poly_evaluate_vec_fast_precomp(
                out, polynomial.raw()->coeffs, polynomial.raw()->length,
                tree_, length_, precision);
        return true;
    }

    bool matches(slong length, slong precision) const noexcept {
        return tree_ != nullptr && length_ == length &&
               precision_ == precision;
    }

    slong length() const noexcept { return length_; }
    slong precision() const noexcept { return precision_; }

private:
    acb_ptr* tree_ = nullptr;
    slong length_ = 0;
    slong precision_ = 0;
};

inline void swap(AcbPoly& left, AcbPoly& right) noexcept { left.swap(right); }

inline void swap(AcbPolyEvaluationTree& left,
                 AcbPolyEvaluationTree& right) noexcept {
    left.swap(right);
}

}  // namespace silex::flint
