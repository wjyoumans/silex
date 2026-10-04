// Tests for the proven generation bound of factor_base_class_group_bound in
// degree n >= 3: the minimum of Zimmert's bound (Zimmert 1981, Satz 2, with
// alpha from Bemerkung 1), Minkowski's bound with (4/pi)^r2, and the former
// exact integer Minkowski form with 2^r2.
//
// The pinned values were computed independently of Silex: Satz 2 evaluated
// with mpmath at 60 digits, and the field discriminants and signatures with
// an independent computer algebra system.

#include <silex/flint/fmpz.hpp>

#include "factor_base/factor_base_internal.hpp"
#include "test_support.hpp"

#include <flint/fmpz.h>

#include <cassert>

namespace {
namespace sflint = silex::flint;

bool set_decimal(sflint::Fmpz& out, const char* decimal) noexcept {
    return fmpz_set_str(out.raw(), decimal, 10) == 0;
}

// Zimmert's bound alone at |d| = 10^30.
struct ZimmertRow {
    slong r1;
    slong r2;
    const char* bound;
};

constexpr ZimmertRow kZimmertRows[] = {
    {3, 0, "215662196762693"},
    {1, 1, "297995988865044"},
    {4, 0, "69179916168941"},
    {2, 1, "102572524138711"},
    {0, 2, "147219529777701"},
    {6, 0, "5315859884779"},
    {0, 3, "21390761573066"},
    {2, 2, "13768462196552"},
    {8, 0, "323813917622"},
    {0, 4, "2593843008062"},
    {10, 0, "17081788300"},
    {5, 4, "2022011626"},
    {0, 10, "1743322"},
    {16, 0, "1509650"},
    {19, 0, "11839"},
    {20, 0, "2309"},
    {0, 9, "20716959"},
    {1, 9, "4405574"},};

// Every class/unit fixture field of degree >= 3 in test/data and a corpus of
// benchmark fields of degrees 3 to 81: absolute discriminant of the maximal
// order, signature, the former bound, and the current bound.
struct FieldRow {
    const char* name;
    slong r1;
    slong r2;
    const char* abs_discriminant;
    const char* former_bound;
    const char* bound;
};

constexpr FieldRow kFieldRows[] = {
    {"cubic_disc81_proven", 3, 0,
     "81",
     "2", "2"},
    {"quartic_x4_plus_1_proven", 0, 2,
     "256",
     "6", "3"},
    {"cubic_disc643_nontrivial_proven", 1, 1,
     "643",
     "12", "8"},
    {"cubic_seeded_h4_proven", 1, 1,
     "744",
     "13", "8"},
    {"quartic_disc892_proven", 0, 2,
     "892",
     "12", "5"},
    {"quartic_disc70640_proven", 2, 1,
     "70640",
     "50", "28"},
    {"quartic_disc223479_proven", 2, 1,
     "223479",
     "89", "49"},
    {"quartic_disc35019_proven", 2, 1,
     "3891",
     "12", "7"},
    {"quartic_disc1412343_proven", 2, 1,
     "1412343",
     "223", "122"},
    {"quintic_disc401370255_proven", 3, 1,
     "401370255",
     "1539", "624"},
    {"cubic_disc23_analytic_grh", 1, 1,
     "23",
     "3", "2"},
    {"cubic_x3_minus_2_proven", 1, 1,
     "108",
     "5", "3"},
    {"cubic_disc2213_proven", 3, 0,
     "2213",
     "11", "11"},
    {"quartic_x4_minus_x_minus_1_proven", 2, 1,
     "283",
     "4", "2"},
    {"quartic_disc1856_proven", 2, 1,
     "1856",
     "9", "5"},
    {"quintic_x5_minus_x_minus_1_proven", 1, 2,
     "2869",
     "9", "3"},
    {"quintic_disc4417_proven", 1, 2,
     "4417",
     "11", "4"},
    {"quintic_disc11119_proven", 3, 1,
     "11119",
     "9", "4"},
    {"quintic_disc57895_proven", 3, 1,
     "57895",
     "19", "8"},
    {"sextic_x6_minus_x_minus_1_proven", 2, 2,
     "49781",
     "14", "4"},
    {"septic_x7_minus_x_minus_1_proven", 1, 3,
     "776887",
     "44", "6"},
    {"octic_x8_minus_x_minus_1_proven", 2, 3,
     "17600759",
     "81", "7"},
    {"nonic_x9_minus_x_minus_1_proven", 1, 4,
     "370643273",
     "289", "13"},
    {"decic_x10_minus_x_minus_1_proven", 2, 4,
     "10387420489",
     "592", "17"},
    {"undecic_x11_minus_x_minus_1_proven", 1, 5,
     "275311670611",
     "2350", "36"},
    {"duodecic_x12_minus_x_minus_1_proven", 2, 5,
     "9201412118867",
     "5215", "49"},
    {"tridecic_x13_minus_x_minus_1_proven", 1, 6,
     "293959006143997",
     "22561", "112"},
    {"tetradecic_x14_minus_x_minus_1_proven", 2, 6,
     "11414881932150269",
     "53646", "158"},
    {"quindecic_x15_minus_x_minus_1_proven", 1, 7,
     "426781883555301359",
     "249715", "389"},
    {"hexadecic_x16_minus_x_minus_1_proven", 2, 7,
     "18884637964090410991",
     "630906", "567"},
    {"septendecic_x17_minus_x_minus_1_proven", 1, 8,
     "808793517812627212561",
     "3130370", "1480"},
    {"octodecic_x18_minus_x_plus_1_proven", 0, 9,
     "38519167813410200811247",
     "16351001", "4066"},
    {"nonadecic_x19_minus_x_minus_1_proven", 1, 9,
     "1939073247585017051548555",
     "43837247", "6135"},
    {"biquadratic_v4_disc2304_proven", 4, 0,
     "2304",
     "5", "4"},
    {"biquadratic_v4_disc7056_proven", 4, 0,
     "7056",
     "8", "6"},
    {"random_d3_h16_s42_proven", 1, 1,
     "155659",
     "176", "112"},
    {"random_d3_h16_s48_proven", 3, 0,
     "97637",
     "70", "68"},
    {"random_d4_h16_s52_proven", 2, 1,
     "64028500",
     "1501", "821"},
    {"random_d4_h16_s54_proven", 0, 2,
     "5251529",
     "860", "338"},
    {"random_d4_h16_s59_proven", 0, 2,
     "465808",
     "257", "101"},
    {"random_d4_h8_s41_proven", 0, 2,
     "65016",
     "96", "38"},
    {"random_d5_h16_s58_proven", 1, 2,
     "432320249",
     "3194", "985"},
    {"random_d5_h16_s43_proven", 3, 1,
     "1744918375",
     "3209", "1301"},
    {"random_d7_h4_s55_proven", 3, 2,
     "100059217373",
     "7744", "1176"},
    {"random_d7_h4_s52_proven", 1, 3,
     "9574084639",
     "4791", "588"},
    {"random_d8_h4_s49_proven", 0, 4,
     "376623820609",
     "23598", "1592"},
    {"random_d8_h2_s48_proven", 2, 3,
     "23168842868",
     "2927", "241"},
    {"random_d6_h4_s11_proven", 2, 2,
     "278655197",
     "1031", "230"},
    {"random_d6_h8_s60_proven", 2, 2,
     "20129259201",
     "8758", "1954"},
    {"random_d6_h8_s64_proven", 2, 2,
     "1096336816",
     "2044", "456"},
    {"random_d6_h4_s12_proven", 4, 1,
     "2631408447",
     "1584", "444"},
    {"random_d9_h4_s8_proven", 3, 3,
     "249476303566836",
     "118355", "6264"},
    {"random_d9_h2_s8_proven", 3, 3,
     "12901951910388",
     "26916", "1425"},
    {"random_d10_h2_s21_proven", 2, 4,
     "251437415419148",
     "92066", "2615"},
    {"random_d10_h4_s27_proven", 0, 5,
     "10072135451714567",
     "1165397", "28185"},
    {"random_d10_h2_s34_proven", 2, 4,
     "12908260015424",
     "20861", "593"},
    {"random_d12_h2_s16_proven", 0, 6,
     "47194289147436245",
     "746942", "6104"},
    {"random_d12_h2_s11_proven", 4, 4,
     "26085893348472064",
     "138831", "1457"},
    {"cyclotomic_5_phi4", 0, 2,
     "125",
     "5", "2"},
    {"cyclotomic_7_phi6", 0, 3,
     "16807",
     "17", "3"},
    {"cyclotomic_9_phi6", 0, 3,
     "19683",
     "18", "4"},
    {"cyclotomic_11_phi10", 0, 5,
     "2357947691",
     "564", "14"},
    {"cyclotomic_12_phi4", 0, 2,
     "144",
     "5", "2"},
    {"cyclotomic_13_phi12", 0, 6,
     "1792160394037",
     "4603", "38"},
    {"cyclotomic_15_phi8", 0, 4,
     "1265625",
     "44", "3"},
    {"cyclotomic_16_phi8", 0, 4,
     "16777216",
     "158", "11"},
    {"cyclotomic_17_phi16", 0, 8,
     "2862423051509815793",
     "491255", "404"},
    {"cyclotomic_19_phi18", 0, 9,
     "5480386857784802185939",
     "6167534", "1534"},
    {"cyclotomic_20_phi8", 0, 4,
     "4000000",
     "77", "6"},
    {"cyclotomic_21_phi12", 0, 6,
     "205924456521",
     "1561", "13"},
    {"cyclotomic_23_phi22", 0, 11,
     "39471584120695485887249589623",
     "1339491139", "9324407"},
    {"cyclotomic_24_phi8", 0, 4,
     "5308416",
     "89", "6"},
    {"cyclotomic_25_phi20", 0, 10,
     "2910383045673370361328125",
     "40532159", "2975"},
    {"cyclotomic_27_phi18", 0, 9,
     "2954312706550833698643",
     "4528289", "1127"},
    {"cyclotomic_28_phi12", 0, 6,
     "1157018619904",
     "3699", "31"},
    {"cyclotomic_32_phi16", 0, 8,
     "18446744073709551616",
     "1247096", "1026"},
    {"cyclotomic_33_phi20", 0, 10,
     "328307557444402776721569",
     "13613353", "999"},
    {"cyclotomic_35_phi24", 0, 12,
     "304383340063522342681884765625",
     "1051250599", "4658733"},
    {"cyclotomic_36_phi12", 0, 6,
     "1586874322944",
     "4332", "36"},
    {"cyclotomic_39_phi24", 0, 12,
     "1706902865139206151939937338729",
     "2489431731", "11032190"},
    {"cyclotomic_40_phi16", 0, 8,
     "1048576000000000000",
     "297331", "245"},
    {"cyclotomic_44_phi20", 0, 10,
     "5829995856912430117421056",
     "57366557", "4210"},
    {"cyclotomic_45_phi24", 0, 12,
     "572565594852444156646728515625",
     "1441811354", "6389545"},
    {"cyclotomic_48_phi16", 0, 8,
     "1846757322198614016",
     "394589", "325"},
    {"cyclotomic_52_phi24", 0, 12,
     "53885714612646242347927893704704",
     "13987259763", "61986072"},
    {"cyclotomic_56_phi24", 0, 12,
     "22459526297810799636782730182656",
     "9030173885", "40018204"},
    {"cyclotomic_60_phi16", 0, 8,
     "104976000000000000",
     "94078", "78"},
    {"cyclotomic_72_phi24", 0, 12,
     "42247883974617233597120303333376",
     "12385065221", "54885772"},
    {"cyclotomic_84_phi24", 0, 12,
     "711435861303500483618465120256",
     "1607176944", "7122381"},
    {"cyclotomic_51_phi32", 0, 16,
     "352701833122210710593389803720131611763844129",
     "221593098346526", "161301422505"},
    {"cyclotomic_64_phi32", 0, 16,
     "1461501637330902918203684832716283019655932542976",
     "14264351252571977", "10383266289959"},
    {"cyclotomic_37_phi36", 0, 18,
     "7710105884424969623139759010953858981831553019262380893",
     "2545160735010903888", "750857335698544"},
    {"cyclotomic_41_phi40", 0, 20,
     "791717805254439023624865699561776475898803884688668051353443"
         "161",
     "1991271758539031907701", "238085525850535145"},
    {"cyclotomic_real_23_deg11", 11, 0,
     "41426511213649",
     "901", "25"},
    {"cyclotomic_real_29_deg14", 14, 0,
     "10260628712958602189",
     "25131", "116"},
    {"cyclotomic_real_31_deg15", 15, 0,
     "756943935220796320321",
     "82161", "205"},
    {"cyclotomic_real_37_deg18", 18, 0,
     "456487940826035155404146917",
     "3476567", "1286"},
    {"cyclotomic_real_41_deg20", 20, 0,
     "4394336169668803158610484050361",
     "48637512", "4839"},
    {"cyclotomic_real_64_deg16", 16, 0,
     "604462909807314587353088",
     "881830", "1174"},
    {"cyclotomic_real_136_deg32", 32, 0,
     "151142765320815856472639544599622796222554813389533609984",
     "2213431686212152", "2213431686212152"},
    {"cyclotomic_real_163_deg81", 81, 0,
     "944079032563526398867950075336198896160191313866593600362194"
         "402467918886055306297936965237879312325253718038737301375553"
         "113182090904662045425093257516971382960144678736277019201",
     "4607129166135942287226553243039933095298503353984207219", "4607129166135942287226553243039933095298503353984207219"},
    {"pure_x3_minus_6", 1, 1,
     "972",
     "15", "9"},
    {"pure_x3_minus_7", 1, 1,
     "1323",
     "17", "11"},
    {"pure_x3_minus_10", 1, 1,
     "300",
     "8", "5"},
    {"pure_x4_minus_2", 2, 1,
     "2048",
     "9", "5"},
    {"pure_x4_minus_3", 2, 1,
     "6912",
     "16", "9"},
    {"pure_x4_minus_5", 2, 1,
     "2000",
     "9", "5"},
    {"pure_x4_minus_7", 2, 1,
     "87808",
     "56", "31"},
    {"pure_x5_minus_2", 1, 2,
     "50000",
     "35", "11"},
    {"pure_x5_minus_3", 1, 2,
     "253125",
     "78", "24"},
    {"pure_x5_minus_7", 1, 2,
     "300125",
     "85", "26"},
    {"pure_x6_minus_2", 2, 2,
     "1492992",
     "76", "17"},
    {"pure_x6_minus_3", 2, 2,
     "11337408",
     "208", "47"},
    {"pure_x6_minus_5", 2, 2,
     "2278125",
     "94", "21"},
    {"pure_x7_minus_2", 1, 3,
     "52706752",
     "356", "44"},
    {"pure_x7_minus_3", 1, 3,
     "600362847",
     "1200", "148"},
    {"pure_x8_minus_2", 2, 3,
     "2147483648",
     "891", "74"},
    {"pure_x8_minus_3", 2, 3,
     "36691771392",
     "3683", "304"},
    {"pure_x8_minus_5", 2, 3,
     "5120000000",
     "1376", "114"},
    {"simplest_cubic_a1", 3, 0,
     "169",
     "3", "3"},
    {"simplest_cubic_a2", 3, 0,
     "361",
     "5", "5"},
    {"simplest_cubic_a4", 3, 0,
     "1369",
     "9", "8"},
    {"simplest_cubic_a5", 3, 0,
     "49",
     "2", "2"},
    {"simplest_cubic_a6", 3, 0,
     "3969",
     "14", "14"},
    {"simplest_cubic_a8", 3, 0,
     "9409",
     "22", "21"},
    {"simplest_cubic_a10", 3, 0,
     "19321",
     "31", "30"},
    {"simplest_quartic_t1", 4, 0,
     "4913",
     "7", "5"},
    {"simplest_quartic_t2", 4, 0,
     "2000",
     "5", "4"},
    {"simplest_quartic_t4", 4, 0,
     "2048",
     "5", "4"},
    {"simplest_quartic_t5", 4, 0,
     "68921",
     "25", "19"},
    {"simplest_quartic_t6", 4, 0,
     "35152",
     "18", "13"},
    {"simplest_quartic_t7", 4, 0,
     "274625",
     "50", "37"},
    {"simplest_quartic_t8", 4, 0,
     "8000",
     "9", "7"},
    {"simplest_quartic_t10", 4, 0,
     "390224",
     "59", "44"},
    {"simplest_quartic_t12", 4, 0,
     "256000",
     "48", "36"},
    {"multiquadratic_m5_p13", 0, 2,
     "67600",
     "98", "39"},
    {"multiquadratic_m23_m31", 0, 2,
     "508369",
     "268", "105"},
    {"multiquadratic_p10_m47", 0, 2,
     "3534400",
     "705", "277"},
    {"multiquadratic_m5_m23", 0, 2,
     "211600",
     "173", "68"},
    {"multiquadratic_m1_m3_p5", 0, 4,
     "12960000",
     "139", "10"},
    {"multiquadratic_m3_m5_m7", 0, 4,
     "31116960000",
     "6783", "458"},
    {"multiquadratic_m23_m31_m47", 0, 4,
     "1261100073931868641",
     "43181268", "2912853"},
    {"multiquadratic_p5_m7_m11", 0, 4,
     "21970650625",
     "5700", "385"},
    {"lmfdb_6_6_46411625_1", 6, 0,
     "46411625",
     "106", "37"},
    {"lmfdb_6_6_49744125_1", 6, 0,
     "49744125",
     "109", "38"},
    {"lmfdb_8_8_73116160000_2", 8, 0,
     "73116160000",
     "650", "88"},
    {"lmfdb_8_8_82944000000_1", 8, 0,
     "82944000000",
     "693", "94"},
    {"lmfdb_8_0_64000000_3", 0, 4,
     "64000000",
     "308", "21"},
    {"lmfdb_8_0_76488065_1", 0, 4,
     "76488065",
     "337", "23"},};

int test_not_above_former_bound() {
    const char* discriminants[] = {
            "1",
            "3",
            "23",
            "1000003",
            "1000000000000",
            "10000000000000000000000000",
            "10000000000000000000000000000000000000000",
            "1606938044258990275541962092341162602522202993782792835301377",
    };
    sflint::Fmpz abs_discriminant;
    sflint::Fmpz bound;
    sflint::Fmpz former;
    for (slong degree = 3; degree <= 24; ++degree) {
        for (slong r2 = 0; 2 * r2 <= degree; ++r2) {
            const slong r1 = degree - 2 * r2;
            for (const char* decimal : discriminants) {
                if (!set_decimal(abs_discriminant, decimal) ||
                    !silex::detail::generation_bound(
                            sflint::FmpzRef(bound),
                            sflint::FmpzConstRef(abs_discriminant), r1, r2)) {
                    return 1;
                }
                silex::test::former_generic_generation_bound(
                        former, abs_discriminant, degree, r2);
                if (fmpz_cmp_ui(bound.raw(), 1) < 0 ||
                    fmpz_cmp(bound.raw(), former.raw()) > 0) {
                    return 1;
                }
            }
            // For large |d| the bound is strictly smaller whenever a
            // sharper constant applies: Zimmert's bound up to degree 20,
            // and (4/pi)^r2 < 2^r2 once r2 > 0.
            if (degree <= 20 || r2 > 0) {
                if (!set_decimal(abs_discriminant, discriminants[6]) ||
                    !silex::detail::generation_bound(
                            sflint::FmpzRef(bound),
                            sflint::FmpzConstRef(abs_discriminant), r1, r2)) {
                    return 1;
                }
                silex::test::former_generic_generation_bound(
                        former, abs_discriminant, degree, r2);
                if (fmpz_cmp(bound.raw(), former.raw()) >= 0) {
                    return 1;
                }
            }
        }
    }
    return 0;
}

int test_zimmert_pinned_values() {
    sflint::Fmpz abs_discriminant;
    sflint::Fmpz bound;
    sflint::Fmpz expected;
    fmpz_set_ui(abs_discriminant.raw(), 10);
    fmpz_pow_ui(abs_discriminant.raw(), abs_discriminant.raw(), 30);
    for (const ZimmertRow& row : kZimmertRows) {
        if (!set_decimal(expected, row.bound) ||
            !silex::detail::zimmert_generation_bound(
                    sflint::FmpzRef(bound),
                    sflint::FmpzConstRef(abs_discriminant), row.r1, row.r2) ||
            !fmpz_equal(bound.raw(), expected.raw())) {
            return 1;
        }
    }
    return 0;
}

int test_field_pinned_values() {
    sflint::Fmpz abs_discriminant;
    sflint::Fmpz bound;
    sflint::Fmpz former;
    sflint::Fmpz expected;
    for (const FieldRow& row : kFieldRows) {
        const slong degree = row.r1 + 2 * row.r2;
        if (!set_decimal(abs_discriminant, row.abs_discriminant) ||
            !silex::detail::generation_bound(
                    sflint::FmpzRef(bound),
                    sflint::FmpzConstRef(abs_discriminant), row.r1, row.r2) ||
            !set_decimal(expected, row.bound) ||
            !fmpz_equal(bound.raw(), expected.raw())) {
            return 1;
        }
        silex::test::former_generic_generation_bound(
                former, abs_discriminant, degree, row.r2);
        if (!set_decimal(expected, row.former_bound) ||
            !fmpz_equal(former.raw(), expected.raw()) ||
            fmpz_cmp(bound.raw(), former.raw()) > 0) {
            return 1;
        }
    }
    return 0;
}

int test_rejected_inputs_preserve_output() {
    sflint::Fmpz abs_discriminant;
    sflint::Fmpz bound;
    fmpz_set_ui(abs_discriminant.raw(), 1000);
    fmpz_set_si(bound.raw(), 99);
    const slong rejected[][2] = {{2, 0}, {0, 1}, {1, 0}, {-1, 2}, {3, -1}};
    for (const auto& signature : rejected) {
        if (silex::detail::generation_bound(
                    sflint::FmpzRef(bound),
                    sflint::FmpzConstRef(abs_discriminant), signature[0],
                    signature[1]) ||
            silex::detail::zimmert_generation_bound(
                    sflint::FmpzRef(bound),
                    sflint::FmpzConstRef(abs_discriminant), signature[0],
                    signature[1])) {
            return 1;
        }
    }
    // Zimmert's bound is tabulated only up to degree 20.
    if (silex::detail::zimmert_generation_bound(
                sflint::FmpzRef(bound),
                sflint::FmpzConstRef(abs_discriminant), 21, 0) ||
        silex::detail::zimmert_generation_bound(
                sflint::FmpzRef(bound),
                sflint::FmpzConstRef(abs_discriminant), 1, 10)) {
        return 1;
    }
    fmpz_zero(abs_discriminant.raw());
    if (silex::detail::generation_bound(
                sflint::FmpzRef(bound),
                sflint::FmpzConstRef(abs_discriminant), 3, 0)) {
        return 1;
    }
    return fmpz_cmp_si(bound.raw(), 99) == 0 ? 0 : 1;
}

}  // namespace

int main() {
    assert(test_not_above_former_bound() == 0);
    assert(test_zimmert_pinned_values() == 0);
    assert(test_field_pinned_values() == 0);
    assert(test_rejected_inputs_preserve_output() == 0);
    return 0;
}
