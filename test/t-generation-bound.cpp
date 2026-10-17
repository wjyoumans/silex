// Tests for the proven generation bound of factor_base_class_group_bound in
// degree n >= 3: the minimum of Zimmert's bound (Zimmert 1981, Satz 2, with
// alpha from Bemerkung 1), Minkowski's bound with (4/pi)^r2, and the former
// exact integer Minkowski form with 2^r2.
//
// The pinned values were computed independently of Silex: Satz 2 evaluated
// with mpmath at 60 digits (degree above 20: at 80 digits by
// test/zimmert_satz2_oracle.py), and the field discriminants and signatures
// with an independent computer algebra system.

#include <silex/factor_base.hpp>
#include <silex/flint/arb.hpp>
#include <silex/flint/arf.hpp>
#include <silex/flint/fmpq_mat.hpp>
#include <silex/flint/fmpq_poly.hpp>
#include <silex/flint/fmpz.hpp>

#include "factor_base/factor_base_internal.hpp"
#include "test_support.hpp"

#include <flint/arb.h>
#include <flint/arf.h>
#include <flint/fmpq_mat.h>
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

// Zimmert's bound alone in degree above 20 at |d| = 10^e, with the gamma
// (hundredths) it uses: test/zimmert_satz2_oracle.py --bounds (mpmath, 80
// digits; gamma maximises R over k/100, 1 <= k <= 300).  Degrees 21 to 64
// take gamma from the generated table, 65 and above from the search.
struct ZimmertHighRow {
    slong r1;
    slong r2;
    ulong exponent;
    ulong gamma_hundredths;
    const char* bound;
};

constexpr ZimmertHighRow kZimmertHighRows[] = {
    {21, 0, 56, 34, "4461816526040420"},
    {1, 10, 50, 44, "3632555759654437"},
    {0, 11, 50, 44, "1428351254775753"},
    {0, 12, 52, 41, "1144949698937169"},
    {32, 0, 72, 27, "4145294854867151"},
    {0, 16, 62, 35, "3954732912154843"},
    {0, 20, 70, 31, "1113085869162007"},
    {64, 0, 120, 18, "1838736847234003"},
    {0, 32, 98, 24, "1185359534516218"},
    {65, 0, 122, 18, "3086360646852981"},
    {1, 32, 100, 23, "2080786071392452"},
    {81, 0, 146, 16, "1047561253645817"},
    {100, 0, 176, 14, "1183658781616385"},
    {0, 50, 142, 18, "4135658928221791"},
};

// Zimmert 1981, Tabelle 3 (p. 377), every row with 3 <= n <= 20: the printed
// gamma (hundredths) and the printed value digits * 10^exponent of
// exp((F + n log pi)/2), rounded down.  Silex's enclosure of that constant
// must lie in [printed, printed + 10^exponent).  The same rows, with the
// margins, are printed by test/zimmert_twin_class_oracle.py.
struct TwinTableRow {
    slong r1;
    slong r2;
    ulong gamma_hundredths;
    ulong digits;
    slong exponent;
};

constexpr TwinTableRow kZimmertTable3Rows[] = {
    {3, 0, 153, 556, -2},  {1, 1, 182, 395, -2},  {4, 0, 115, 185, -1},
    {2, 1, 131, 121, -1},  {0, 2, 153, 832, -2},  {5, 0, 95, 685, -1},
    {3, 1, 106, 427, -1},  {1, 2, 119, 274, -1},  {6, 0, 82, 273, 0},
    {0, 3, 109, 629, -1},  {8, 0, 66, 511, 1},    {0, 4, 87, 570, 0},
    {10, 0, 56, 110, 3},   {0, 5, 74, 577, 1},    {20, 0, 36, 163, 10},
    {0, 10, 47, 148, 7},
};

// The twin-class generation bound at |d| = 10^30, from
// test/zimmert_twin_class_oracle.py --bounds (mpmath, 40 digits).
constexpr ZimmertRow kZimmertTwinRows[] = {
    {3, 0, "179543323431641"},
    {1, 1, "253040414206487"},
    {4, 0, "54041713443025"},
    {2, 1, "82010497502656"},
    {0, 2, "120099755237859"},
    {6, 0, "3653623151319"},
    {0, 3, "15874536942412"},
    {8, 0, "195442517658"},
    {0, 4, "1753661445481"},
    {10, 0, "9031932192"},
    {0, 10, "671356"},
    {16, 0, "529770"},
    {20, 0, "611"},
    {1, 9, "1750317"},
};

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
     "1339491139", "28378"},
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
     "1051250599", "6317"},
    {"cyclotomic_36_phi12", 0, 6,
     "1586874322944",
     "4332", "36"},
    {"cyclotomic_39_phi24", 0, 12,
     "1706902865139206151939937338729",
     "2489431731", "14959"},
    {"cyclotomic_40_phi16", 0, 8,
     "1048576000000000000",
     "297331", "245"},
    {"cyclotomic_44_phi20", 0, 10,
     "5829995856912430117421056",
     "57366557", "4210"},
    {"cyclotomic_45_phi24", 0, 12,
     "572565594852444156646728515625",
     "1441811354", "8664"},
    {"cyclotomic_48_phi16", 0, 8,
     "1846757322198614016",
     "394589", "325"},
    {"cyclotomic_52_phi24", 0, 12,
     "53885714612646242347927893704704",
     "13987259763", "84048"},
    {"cyclotomic_56_phi24", 0, 12,
     "22459526297810799636782730182656",
     "9030173885", "54261"},
    {"cyclotomic_60_phi16", 0, 8,
     "104976000000000000",
     "94078", "78"},
    {"cyclotomic_72_phi24", 0, 12,
     "42247883974617233597120303333376",
     "12385065221", "74420"},
    {"cyclotomic_84_phi24", 0, 12,
     "711435861303500483618465120256",
     "1607176944", "9658"},
    {"cyclotomic_51_phi32", 0, 16,
     "352701833122210710593389803720131611763844129",
     "221593098346526", "7427130"},
    {"cyclotomic_64_phi32", 0, 16,
     "1461501637330902918203684832716283019655932542976",
     "14264351252571977", "478097873"},
    {"cyclotomic_37_phi36", 0, 18,
     "7710105884424969623139759010953858981831553019262380893",
     "2545160735010903888", "5952820917"},
    {"cyclotomic_41_phi40", 0, 20,
     "791717805254439023624865699561776475898803884688668051353443"
         "161",
     "1991271758539031907701", "313194318666"},
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
     "2213431686212152", "50962311"},
    {"cyclotomic_real_163_deg81", 81, 0,
     "944079032563526398867950075336198896160191313866593600362194"
         "402467918886055306297936965237879312325253718038737301375553"
         "113182090904662045425093257516971382960144678736277019201",
     "4607129166135942287226553243039933095298503353984207219", "3218722996420627348246212878036"},
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
            // For large |d| the bound is strictly smaller: Zimmert's bound
            // applies in every degree.
            {
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

// Scales value by 10^(-exponent), exactly for exponent <= 0.
void scale_by_power_of_ten(sflint::Arb& value, slong exponent) noexcept {
    sflint::Arb power;
    arb_set_ui(power.raw(), 10);
    arb_pow_ui(power.raw(), power.raw(),
               static_cast<ulong>(exponent < 0 ? -exponent : exponent), 256);
    if (exponent < 0) {
        arb_mul(value.raw(), value.raw(), power.raw(), 256);
    } else {
        arb_div(value.raw(), value.raw(), power.raw(), 256);
    }
}

int test_zimmert_twin_class_table3() {
    sflint::Arb constant;
    sflint::Arb value;
    sflint::Arb edge;
    for (const TwinTableRow& row : kZimmertTable3Rows) {
        if (!silex::detail::zimmert_twin_class_log_lower_bound(
                    sflint::ArbRef(constant), row.r1, row.r2,
                    row.gamma_hundredths)) {
            return 1;
        }
        arb_exp(value.raw(), constant.raw(), 256);
        scale_by_power_of_ten(value, row.exponent);
        arb_set_ui(edge.raw(), row.digits);
        if (arb_ge(value.raw(), edge.raw()) == 0) {
            return 1;
        }
        arb_add_ui(edge.raw(), edge.raw(), 1, 256);
        if (arb_lt(value.raw(), edge.raw()) == 0) {
            return 1;
        }
    }
    return 0;
}

int test_zimmert_high_degree_pinned_values() {
    sflint::Fmpz abs_discriminant;
    sflint::Fmpz bound;
    sflint::Fmpz combined;
    sflint::Fmpz expected;
    for (const ZimmertHighRow& row : kZimmertHighRows) {
        fmpz_ui_pow_ui(abs_discriminant.raw(), 10, row.exponent);
        if (silex::detail::zimmert_satz2_gamma_hundredths(row.r1, row.r2) !=
                    row.gamma_hundredths ||
            !set_decimal(expected, row.bound) ||
            !silex::detail::zimmert_generation_bound(
                    sflint::FmpzRef(bound),
                    sflint::FmpzConstRef(abs_discriminant), row.r1, row.r2) ||
            !fmpz_equal(bound.raw(), expected.raw())) {
            return 1;
        }
        // Above degree 20 Satz 2 is far below both Minkowski forms, so it
        // is the combined bound.
        if (!silex::detail::generation_bound(
                    sflint::FmpzRef(combined),
                    sflint::FmpzConstRef(abs_discriminant), row.r1, row.r2) ||
            !fmpz_equal(combined.raw(), expected.raw())) {
            return 1;
        }
    }
    return 0;
}

// Zimmert 1981, Tabelle 1 (p. 368), the rows with n = 20 and n = 100: the
// printed gamma (hundredths) and the printed value digits * 10^exponent of
// d^(1/2) / N(a) = exp(R), rounded down to four digits.  Silex's enclosure
// of exp(R) at that gamma must lie in [printed, printed + 10^exponent), and
// Silex's choice of gamma must be the printed one.
constexpr TwinTableRow kZimmertTable1Rows[] = {
    {20, 0, 36, 4332, 8},
    {0, 10, 46, 5736, 5},
    {100, 0, 14, 8448, 69},
    {0, 50, 18, 2417, 52},
};

int test_zimmert_table1() {
    sflint::Arb constant;
    sflint::Arb value;
    sflint::Arb edge;
    for (const TwinTableRow& row : kZimmertTable1Rows) {
        if (silex::detail::zimmert_satz2_gamma_hundredths(row.r1, row.r2) !=
                    row.gamma_hundredths ||
            !silex::detail::zimmert_satz2_log_lower_bound(
                    sflint::ArbRef(constant), row.r1, row.r2,
                    row.gamma_hundredths)) {
            return 1;
        }
        arb_exp(value.raw(), constant.raw(), 256);
        scale_by_power_of_ten(value, row.exponent);
        arb_set_ui(edge.raw(), row.digits);
        if (arb_ge(value.raw(), edge.raw()) == 0) {
            return 1;
        }
        arb_add_ui(edge.raw(), edge.raw(), 1, 256);
        if (arb_lt(value.raw(), edge.raw()) == 0) {
            return 1;
        }
    }
    return 0;
}

// Midpoint of the Satz 2 enclosure at gamma = k/100; false if k is outside
// the search grid or the enclosure fails.
bool satz2_midpoint(sflint::Arf& out, slong r1, slong r2, ulong k) noexcept {
    sflint::Arb value;
    if (k == 0 || k > 300 ||
        !silex::detail::zimmert_satz2_log_lower_bound(sflint::ArbRef(value),
                                                      r1, r2, k)) {
        return false;
    }
    arf_set(out.raw(), arb_midref(value.raw()));
    return true;
}

// The gamma of Satz 2: the hand-checked table (n <= 20) and the generated
// table (21 <= n <= 64) agree with the search, and the search is used above.
// The search is run in full for every n <= 20 and for sample degrees in the
// generated range; every other generated entry k is checked to be the grid
// argmax locally (midpoint at k above that at k - 1 and not below that at
// k + 1).  The whole generated table also equals the independent mpmath
// argmax of test/zimmert_satz2_oracle.py --gammas 64 --check.
int test_zimmert_gamma_choice() {
    sflint::Arf here;
    sflint::Arf neighbour;
    for (slong degree = 3; degree <= 66; ++degree) {
        const bool full = degree <= 20 || degree == 21 || degree == 22 ||
                          degree == 43 || degree >= 64;
        for (slong r2 = 0; 2 * r2 <= degree; ++r2) {
            const slong r1 = degree - 2 * r2;
            const ulong gamma =
                    silex::detail::zimmert_satz2_gamma_hundredths(r1, r2);
            if (gamma == 0 || gamma > 300) {
                return 1;
            }
            if (full) {
                if (silex::detail::zimmert_satz2_gamma_search(r1, r2) !=
                    gamma) {
                    return 1;
                }
                continue;
            }
            if (!satz2_midpoint(here, r1, r2, gamma)) {
                return 1;
            }
            if (gamma > 1 &&
                (!satz2_midpoint(neighbour, r1, r2, gamma - 1) ||
                 arf_cmp(here.raw(), neighbour.raw()) <= 0)) {
                return 1;
            }
            if (gamma < 300 &&
                (!satz2_midpoint(neighbour, r1, r2, gamma + 1) ||
                 arf_cmp(here.raw(), neighbour.raw()) < 0)) {
                return 1;
            }
        }
    }
    return 0;
}

int test_zimmert_twin_class_pinned_values() {
    sflint::Fmpz abs_discriminant;
    sflint::Fmpz bound;
    sflint::Fmpz expected;
    fmpz_set_ui(abs_discriminant.raw(), 10);
    fmpz_pow_ui(abs_discriminant.raw(), abs_discriminant.raw(), 30);
    for (const ZimmertRow& row : kZimmertTwinRows) {
        if (!set_decimal(expected, row.bound) ||
            !silex::detail::twin_class_generation_bound(
                    sflint::FmpzRef(bound),
                    sflint::FmpzConstRef(abs_discriminant), row.r1, row.r2) ||
            !fmpz_equal(bound.raw(), expected.raw())) {
            return 1;
        }
    }
    return 0;
}

// The precomputed twin-class table (src/factor_base/zimmert_twin_class_table.inc,
// written by tools/zimmert_twin_class_table.cpp).  For every signature with
// 3 <= n <= 20 the entry exists, carries the tabulated gamma, and its exact
// value E = s / 10^k satisfies U <= E <= U (1 + 10^-kTwinTableTolerance),
// where U is the upper endpoint of the 128-bit Arb enclosure of exp(-c)
// recomputed here from zimmert_twin_class_log_lower_bound, the certified
// enclosure of c = (F + n log pi)/2 (GB-zimmert-F-tail-bound).
// E >= U >= exp(-c) is the safe side: the run-time bound ceil(E sqrt|d|) is
// then at least the theorem's sqrt|d| exp(-c).  The generator rounds U up to
// 30 significant digits (relative gap below 10^-29; the enclosure itself is
// accurate to about 2^-120), so 10^-25 leaves room for a slightly different
// enclosure
// radius under another FLINT build while keeping the table tight:
// ceil(E sqrt|d|) can exceed ceil(U sqrt|d|), the value the former run-time
// Arb evaluation produced up to its own rounding, only if an integer lies in
// [U sqrt|d|, E sqrt|d|], an interval shorter than 10^-5 for every bound
// below 10^20.
constexpr slong kTwinTablePrecision = 128;  // as in the generator
constexpr ulong kTwinTableTolerance = 25;

int test_zimmert_twin_class_verified_table() {
    sflint::Fmpz significand;
    sflint::Fmpz power;
    sflint::Fmpz tolerance;
    sflint::Arb constant;
    sflint::Arb value;
    sflint::Arf upper;
    sflint::Arf scaled;
    sflint::Arf entry;
    sflint::Arf gap;
    fmpz_ui_pow_ui(tolerance.raw(), 10, kTwinTableTolerance);
    slong entries = 0;
    for (slong degree = 2; degree <= 21; ++degree) {
        for (slong r2 = 0; 2 * r2 <= degree; ++r2) {
            const slong r1 = degree - 2 * r2;
            const ulong gamma =
                    silex::detail::zimmert_twin_class_gamma_hundredths(r1, r2);
            slong scale = -1;
            ulong entry_gamma = 0;
            const bool has_entry =
                    silex::detail::zimmert_twin_class_table_entry(
                            sflint::FmpzRef(significand), scale, entry_gamma,
                            r1, r2);
            if (degree < 3 || degree > 20) {
                if (gamma != 0 || has_entry) {
                    return 1;
                }
                continue;
            }
            if (gamma == 0 || !has_entry || entry_gamma != gamma ||
                scale < 0 || fmpz_sgn(significand.raw()) <= 0 ||
                !silex::detail::zimmert_twin_class_log_lower_bound(
                        sflint::ArbRef(constant), r1, r2, gamma)) {
                return 1;
            }
            ++entries;
            arb_neg(constant.raw(), constant.raw());
            arb_exp(value.raw(), constant.raw(), kTwinTablePrecision);
            if (arb_is_finite(value.raw()) == 0) {
                return 1;
            }
            arb_get_ubound_arf(upper.raw(), value.raw(), kTwinTablePrecision);
            // Compare s with U 10^k, both exactly.
            fmpz_ui_pow_ui(power.raw(), 10, static_cast<ulong>(scale));
            arf_mul_fmpz(scaled.raw(), upper.raw(), power.raw(),
                         ARF_PREC_EXACT, ARF_RND_CEIL);
            arf_set_fmpz(entry.raw(), significand.raw());
            if (arf_cmp(entry.raw(), scaled.raw()) < 0) {
                return 1;
            }
            arf_sub(gap.raw(), entry.raw(), scaled.raw(), ARF_PREC_EXACT,
                    ARF_RND_CEIL);
            arf_mul_fmpz(gap.raw(), gap.raw(), tolerance.raw(),
                         ARF_PREC_EXACT, ARF_RND_CEIL);
            if (arf_cmp(gap.raw(), scaled.raw()) > 0) {
                return 1;
            }
        }
    }
    return entries == 117 ? 0 : 1;
}

// The twin-class bound never exceeds generation_bound, is defined exactly
// for 3 <= n <= 20, and rejects gamma outside (0, 3).  Integer gamma is
// accepted: F extends continuously there (GB-zimmert-F-tail-bound, part
// (iii)), and the enclosure needs no special case.
int test_zimmert_twin_class_domain() {
    const char* discriminants[] = {"23", "1000", "1000000007",
                                   "100000000000000000000"};
    sflint::Fmpz abs_discriminant;
    sflint::Fmpz bound;
    sflint::Fmpz classical;
    for (slong degree = 3; degree <= 21; ++degree) {
        for (slong r2 = 0; 2 * r2 <= degree; ++r2) {
            const slong r1 = degree - 2 * r2;
            for (const char* decimal : discriminants) {
                fmpz_set_si(bound.raw(), -7);
                if (!set_decimal(abs_discriminant, decimal) ||
                    !silex::detail::generation_bound(
                            sflint::FmpzRef(classical),
                            sflint::FmpzConstRef(abs_discriminant), r1, r2)) {
                    return 1;
                }
                const bool ok = silex::detail::twin_class_generation_bound(
                        sflint::FmpzRef(bound),
                        sflint::FmpzConstRef(abs_discriminant), r1, r2);
                if (degree > 20) {
                    if (ok || fmpz_cmp_si(bound.raw(), -7) != 0) {
                        return 1;
                    }
                    continue;
                }
                if (!ok || fmpz_cmp_ui(bound.raw(), 1) < 0 ||
                    fmpz_cmp(bound.raw(), classical.raw()) > 0) {
                    return 1;
                }
            }
        }
    }
    sflint::Arb constant;
    const ulong rejected_gammas[] = {0, 300, 450};
    for (ulong gamma : rejected_gammas) {
        if (silex::detail::zimmert_twin_class_log_lower_bound(
                    sflint::ArbRef(constant), 3, 0, gamma)) {
            return 1;
        }
    }
    const ulong integer_gammas[] = {100, 200};
    for (ulong gamma : integer_gammas) {
        if (!silex::detail::zimmert_twin_class_log_lower_bound(
                    sflint::ArbRef(constant), 3, 0, gamma)) {
            return 1;
        }
    }
    if (silex::detail::zimmert_twin_class_log_lower_bound(
                sflint::ArbRef(constant), 0, 0, 50) ||
        !silex::detail::zimmert_twin_class_log_lower_bound(
                sflint::ArbRef(constant), 0, 2, 299)) {
        return 1;
    }
    fmpz_set_ui(abs_discriminant.raw(), 1000);
    if (silex::detail::twin_class_generation_bound(
                sflint::FmpzRef(bound),
                sflint::FmpzConstRef(abs_discriminant), 2, 0)) {
        return 1;
    }
    return 0;
}

// Zimmert's F_{a,b}(gamma) (p. 376) as enclosed by zimmert_f_enclosure
// (GB-zimmert-F-tail-bound, part (iv)).  At non-integer gamma the expected
// values are Zimmert's series summed directly by
// test/zimmert_twin_class_oracle.py --self-check (mpmath nsum, no tail
// estimate; nsum and Levin's transformation agree); at integer gamma, where
// F is the continuous extension, they are the balls of silex-research's
// zimmert_F (src/silex_research/zimmert.py at commit 06191ff, 128 bits).
struct FEnclosureRow {
    slong a;
    slong b;
    slong numerator;
    slong denominator;
    const char* value;
    slong error_bits;  // |value - F| < 2^-error_bits
};

constexpr FEnclosureRow kFEnclosureRows[] = {
    {3, 0, 153, 100, "0.000487830472574483870829948424647818842197", 120},
    {1, 1, 182, 100, "-1.629146269721465997389016536729880835443", 120},
    {20, 0, 36, 100, "33.351959734715381368703448116", 85},
    {2, 2, 230, 100, "-0.586458985898003158518751793902", 95},
    {2, 0, 1, 1, "-1.81257129920006115617569433814368101565", 120},
    {1, 1, 2, 1, "-1.54626400006438097589207069697944921293", 120},
};

int test_zimmert_f_enclosure() {
    sflint::Arb value;
    sflint::Arb expected;
    sflint::Arb coarse;
    for (const FEnclosureRow& row : kFEnclosureRows) {
        slong terms = 0;
        slong order = 0;
        if (!silex::detail::zimmert_f_enclosure(
                    sflint::ArbRef(value), terms, order, row.a, row.b,
                    row.numerator, row.denominator) ||
            terms < 1 || order < 2 ||
            arb_rel_accuracy_bits(value.raw()) < 100) {
            return 1;
        }
        arb_set_str(expected.raw(), row.value, 256);
        arb_add_error_2exp_si(expected.raw(), -row.error_bits);
        if (arb_overlaps(value.raw(), expected.raw()) == 0) {
            return 1;
        }
        // F does not depend on the admissible (L, K) used: a short head
        // and a low order give a wider ball that still contains it.
        terms = 40;
        order = 6;
        if (!silex::detail::zimmert_f_enclosure(
                    sflint::ArbRef(coarse), terms, order, row.a, row.b,
                    row.numerator, row.denominator) ||
            terms != 40 || order != 6 ||
            arb_contains(coarse.raw(), value.raw()) == 0 ||
            arb_rel_accuracy_bits(coarse.raw()) >= 100) {
            return 1;
        }
    }
    // Inadmissible parameters at gamma = 23/10 (c = 28/5, rho = 15/56,
    // sigma = 15/112): L must exceed gamma/2, K must be at least 2, and
    // 2L + 2 must exceed (K - 1)/sigma.
    const slong rejected[][2] = {{1, 6}, {40, 1}, {40, 12}};
    for (const auto& parameters : rejected) {
        slong terms = parameters[0];
        slong order = parameters[1];
        if (silex::detail::zimmert_f_enclosure(sflint::ArbRef(value), terms,
                                               order, 2, 2, 23, 10)) {
            return 1;
        }
    }
    slong terms = 4;  // 2L + 2 = 10 > 1/sigma = 112/15
    slong order = 2;
    if (!silex::detail::zimmert_f_enclosure(sflint::ArbRef(value), terms,
                                            order, 2, 2, 23, 10)) {
        return 1;
    }
    // a >= 1, b >= 0 and gamma > 0 are required.
    terms = 0;
    order = 0;
    if (silex::detail::zimmert_f_enclosure(sflint::ArbRef(value), terms,
                                           order, 0, 1, 1, 2) ||
        silex::detail::zimmert_f_enclosure(sflint::ArbRef(value), terms,
                                           order, 1, -1, 1, 2) ||
        silex::detail::zimmert_f_enclosure(sflint::ArbRef(value), terms,
                                           order, 1, 0, 0, 1)) {
        return 1;
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
    // The gamma of Satz 2 and its enclosure reject the same signatures.
    sflint::Arb log_lower;
    arb_set_si(log_lower.raw(), 7);
    for (const auto& signature : rejected) {
        if (silex::detail::zimmert_satz2_gamma_hundredths(signature[0],
                                                         signature[1]) != 0 ||
            silex::detail::zimmert_satz2_gamma_search(signature[0],
                                                      signature[1]) != 0 ||
            silex::detail::zimmert_satz2_log_lower_bound(
                    sflint::ArbRef(log_lower), signature[0], signature[1],
                    50)) {
            return 1;
        }
    }
    if (silex::detail::zimmert_satz2_log_lower_bound(
                sflint::ArbRef(log_lower), 3, 0, 0) ||
        arb_equal_si(log_lower.raw(), 7) == 0) {
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

// The maximal order of the field defined by
// c3 x^3 + c1 x + c0, computed from the order with basis 1, c3 beta,
// c3^2 beta^2, which is Z[c3 beta] and so lies in O_K.
silex::Order cubic_maximal_order(slong c3, slong c1, slong c0) noexcept {
    sflint::FmpqPoly polynomial;
    sflint::fmpq_poly_zero(polynomial);
    sflint::fmpq_poly_set_coeff_si(polynomial, 3, c3);
    sflint::fmpq_poly_set_coeff_si(polynomial, 1, c1);
    sflint::fmpq_poly_set_coeff_si(polynomial, 0, c0);
    silex::NumberField field = silex::test::field_by_polynomial(
            sflint::FmpqPolyConstRef(polynomial));
    sflint::FmpqMat basis(3, 3);
    fmpq_set_si(fmpq_mat_entry(basis.raw(), 0, 0), 1, 1);
    fmpq_set_si(fmpq_mat_entry(basis.raw(), 1, 1), c3, 1);
    fmpq_set_si(fmpq_mat_entry(basis.raw(), 2, 2), c3 * c3, 1);
    silex::Order integral =
            silex::Order::from_basis(field, sflint::FmpqMatConstRef(basis));
    assert(integral.is_defined());
    silex::Order maximal(field);
    assert(maximal.maximal_order(integral));
    return maximal;
}

// The proven generation requirement on five presentations of the field of
// x^3 + x + 200 (|d| = 1080004, signature (1, 1), class group Z/2; indices
// and prime norms from an independent computer algebra system), where the
// Minkowski-type bound is ceil(6/27 (4/pi) sqrt(1080004)) = 295 and
// Zimmert's twin-class bound is ceil(sqrt(1080004) / 3.95193...) = 263.
// Zimmert's Korollar 1 bound generates the class group with no hypothesis
// on the class of the different (GB-kor1-generation), so every
// presentation gets requirement 263 with route kor1_generation, whatever
// the polynomial order:
// - x^3 + x + 200 has index 1;
// - x^3 + 25 x + 25000 (theta = 5 alpha) has index 125, every prime above
//   5 of norm 5;
// - x^3 + 289 x + 982600 (theta = 17 alpha) has index 17^3, with a prime
//   of norm 289 above 17 (the former conductor_primes route gave 289);
// - x^3 + 361 x + 1371800 (theta = 19 alpha) has a prime of norm 361 above
//   19, above the Minkowski-type bound (formerly no route, 295);
// - 4 x^3 + x + 100 (beta = alpha / 2) is not monic (formerly no route,
//   295).
int test_proven_generation_requirement() {
    struct Row {
        slong c3;
        slong c1;
        slong c0;
        slong requirement;
        silex::detail::GenerationDifferentRoute route;
    };
    using Route = silex::detail::GenerationDifferentRoute;
    const Row rows[] = {
            {1, 1, 200, 263, Route::kor1_generation},
            {1, 25, 25000, 263, Route::kor1_generation},
            {1, 289, 982600, 263, Route::kor1_generation},
            {1, 361, 1371800, 263, Route::kor1_generation},
            {4, 1, 100, 263, Route::kor1_generation},
    };
    sflint::Fmpz bound;
    sflint::Fmpz classical;
    sflint::Fmpz discriminant;
    for (const Row& row : rows) {
        const silex::Order order =
                cubic_maximal_order(row.c3, row.c1, row.c0);
        if (!order.discriminant(sflint::FmpzRef(discriminant)) ||
            !sflint::fmpz_equal_si(discriminant, -1080004) ||
            !silex::factor_base_class_group_bound(sflint::FmpzRef(classical),
                                                  order) ||
            !sflint::fmpz_equal_si(classical, 295)) {
            return 1;
        }
        Route route = Route::none;
        if (!silex::detail::proven_generation_requirement(
                    sflint::FmpzRef(bound), route, order) ||
            !sflint::fmpz_equal_si(bound, row.requirement) ||
            route != row.route) {
            return 1;
        }
    }

    // Degree 2 keeps factor_base_class_group_bound.
    silex::NumberField quadratic = silex::test::quadratic_field(-23);
    silex::Order quadratic_maximal(quadratic);
    if (!quadratic_maximal.maximal_order(
                silex::test::equation_order(quadratic)) ||
        !silex::factor_base_class_group_bound(sflint::FmpzRef(classical),
                                              quadratic_maximal)) {
        return 1;
    }
    Route route = Route::kor1_generation;
    if (!silex::detail::proven_generation_requirement(
                sflint::FmpzRef(bound), route, quadratic_maximal) ||
        !fmpz_equal(bound.raw(), classical.raw()) || route != Route::none) {
        return 1;
    }

    // A non-maximal order has no requirement and leaves the outputs as they
    // were.
    silex::Order equation = silex::test::equation_order(
            *cubic_maximal_order(1, 25, 25000).parent());
    fmpz_set_si(bound.raw(), 99);
    route = Route::kor1_generation;  // a value the function never writes on failure
    if (silex::detail::proven_generation_requirement(sflint::FmpzRef(bound),
                                                     route, equation) ||
        !sflint::fmpz_equal_si(bound, 99) ||
        route != Route::kor1_generation) {
        return 1;
    }
    return 0;
}

}  // namespace

int main() {
    assert(test_not_above_former_bound() == 0);
    assert(test_zimmert_pinned_values() == 0);
    assert(test_zimmert_high_degree_pinned_values() == 0);
    assert(test_zimmert_table1() == 0);
    assert(test_zimmert_gamma_choice() == 0);
    assert(test_field_pinned_values() == 0);
    assert(test_rejected_inputs_preserve_output() == 0);
    assert(test_zimmert_twin_class_table3() == 0);
    assert(test_zimmert_twin_class_pinned_values() == 0);
    assert(test_zimmert_twin_class_verified_table() == 0);
    assert(test_zimmert_twin_class_domain() == 0);
    assert(test_zimmert_f_enclosure() == 0);
    assert(test_proven_generation_requirement() == 0);
    return 0;
}
