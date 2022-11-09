/* When we relax system in new volume we have a problem of PBC of old volume, i.e. we have hydrogens of one water
 * molecule distanced from oxygen to characteristic size of old box. It creates a problem of too rapid moving oxygen and
 * hydrogens to each other. (LAMMPS understands molecules in tip4p.)
 * The solution, that is ee now, is to find the molecules with this problem and replace hydrogens to their oxygen before
 * staring the simulation. For excluding overlap between atoms (overlap creates the same problem) we will generate
 * random angles for offsets of O-H vector until there will be no overlaps. */

#include <iostream>
#include <random>
#include <vector>
#include <fstream>
#include <string>
#include <tuple>
#include <sstream>
#include <cmath>
#include <iterator>


const int first_inserted_atom = 648;
const double OH_dist = 1.0;
const double new_box_size = 19.235;


std::random_device rd;  // Will be used to obtain a seed for the random number engine.
std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd().


typedef std::tuple<double, double, double> data_tuple;


std::vector<data_tuple> coordinates_read (const std::string & name);

void fix (std::vector<data_tuple> & coordinates, const int & first_insert, const double & dist);

void data_file_creation (const std::string & name, std::vector<data_tuple> & data);


int main() {
    std::vector<data_tuple> data = std::move(coordinates_read("coordinates"));
    fix (data, first_inserted_atom, OH_dist);
    data_file_creation("fixed_coordinates", data);
    return 0;
}


namespace std {
    istream& operator >> (istream& in, data_tuple & coordinates) {
        double first, second, third;
        in >> first >> second >> third;
        coordinates = {first, second, third};
        return in;
    }

    ostream& operator << (ostream& out, const data_tuple & coordinates) {
        auto [first, second, third] = coordinates;
        out << first << ' ' << second << ' ' << third << ' ';
        return out;
    }
}


// Read coordinates from text-file.
std::vector<data_tuple> coordinates_read (const std::string & name) {
    std::ifstream fin(name);
    if (!fin.is_open()) throw std::runtime_error("Error opening file.");
    std::vector<data_tuple> tuples_vector;
    copy(std::istream_iterator<data_tuple> {fin},
         std::istream_iterator<data_tuple> {},
         back_inserter(tuples_vector));
    //copy(tuples_vector.begin(), tuples_vector.end(), std::ostream_iterator<data>(std::cout, "\n"));
    return tuples_vector;
}


// Returns distance between two atoms.
template<typename T, size_t... Is>
double distance_impl (T const& t, T const& t1, std::index_sequence<Is...>, std::index_sequence<Is...>) {
    return (std::sqrt((std::pow(std::get<Is>(t) - std::get<Is>(t1), 2) + ...)));
}

template <class Tuple>
double distance (const Tuple & t, const Tuple & t1) {
    constexpr auto size = std::tuple_size<Tuple>{};
    return distance_impl(t, t1, std::make_index_sequence<size>{}, std::make_index_sequence<size>{});
}


// Returns sum of two vectors.
template<size_t Is = 0, typename... Tp>
void vector_offset (std::tuple<Tp...>& vector, std::tuple<Tp...>& offset) {
    std::get<Is>(vector) += std::get<Is>(offset);
    if constexpr(Is + 1 != sizeof...(Tp))
        vector_offset<Is + 1>(vector, offset);
}


// Returns true if atom inside the box.
template<typename T, size_t... Is>
bool out_of_box_impl (T const& t, std::index_sequence<Is...>, const double & box_size) {
    return ((std::get<Is>(t) > box_size) & ...);
}

template <class Tuple>
bool out_of_box (const Tuple & t, const double & box_size) {
    constexpr auto size = std::tuple_size<Tuple>{};
    return out_of_box_impl(t, std::make_index_sequence<size>{}, box_size);
}


// Generate random direction cosines then offsets the atom around ref.
data_tuple rotation (data_tuple coordinate, const double & dist) {
    std::uniform_real_distribution<> dis(0.0, 1.0);
    std::uniform_int_distribution<> dis_int(1, 2);
    double cos_phi, a, b, cos_psi, cos_gamma, k, d = 10; //Why d = 10? Works! Do not touch!
    do {
        do {
            cos_phi = 2 * dis(gen) - 1;
            do {
                a = 2 * dis(gen) - 1;
                b = 2 * dis(gen) - 1;
                d = std::pow(a, 2) + std::pow(b, 2);
            } while (d > 1);
            cos_psi = a / std::sqrt(d);
        } while (std::pow(cos_phi, 2) + std::pow(cos_psi, 2) > 1.0);
        k = dis_int(gen);
        cos_gamma = std::pow(-1, k) * std::sqrt(1.0 - (std::pow(cos_phi, 2) + std::pow(cos_psi, 2)));
        data_tuple offset = std::make_tuple(cos_gamma * dist, cos_phi * dist, cos_psi * dist);
        vector_offset(coordinate, offset);
    } while (!out_of_box(coordinate, new_box_size)); // We have to stay inside the box.
    return coordinate;
}


bool is_equal (const double & a, const double & b) {
    return std::fabs(a - b) < std::numeric_limits<double>::epsilon();
}


template<typename T, size_t... Is>
bool equal_tuples_impl (T const& t, T const& t1, std::index_sequence<Is...>, std::index_sequence<Is...>) {
    return ((is_equal(std::get<Is>(t), std::get<Is>(t1))) & ...);
}

// Returns true if two tuples (t, t1) contains the same numbers.
template <class Tuple>
bool equal_tuples (const Tuple& t, const Tuple& t1) {
    constexpr auto size = std::tuple_size<Tuple>{};
    return equal_tuples_impl(t, t1, std::make_index_sequence<size>{}, std::make_index_sequence<size>{});
}


// Returns true if there's overlap.
bool intersected (data_tuple & atom, std::vector<data_tuple> & coordinates) {
    for (auto & coordinate : coordinates)
        if (distance(atom, coordinate) < OH_dist && !equal_tuples(atom, coordinate))
            return true;
    return false; //Use any_of instead. I just too lazy.
}


void fix (std::vector<data_tuple> & coordinates, const int & first_insert, const double & dist) {
    // NOTE! Distance must be calculated only between H and O!
    for (int i = 0; i < first_insert; i += 3) {
        for (int j = 1; j < 3; ++j) {
            int num = i+j;
            if (distance(coordinates[i], coordinates[num]) > dist) {
                data_tuple buf = coordinates[num];
                while (intersected(buf, coordinates))
                    buf = std::move(rotation(coordinates[num], dist));
                coordinates[num] = buf;
            }
        }
    }
}


// std::to_string not safe enough. It will be used everywhere instead of std::to_string.
template <typename T>
std::string toString (T val) {
    std::ostringstream oss;
    oss << val;
    return oss.str();
}


// Returns string contains tuple content.
template<typename T, size_t... Is>
std::string tuple_to_string_impl (T const& t, std::index_sequence<Is...>) {
    return ((toString(std::get<Is>(t)) + '\t') + ...);
}

template <class Tuple>
std::string tuple_to_string (const Tuple& t) {
    constexpr auto size = std::tuple_size<Tuple>{};
    return tuple_to_string_impl(t, std::make_index_sequence<size>{});
}


// Creates text-file with coordinates from std::vector<data_tuple> with given name.
void data_file_creation (const std::string & name, std::vector<data_tuple> & data) {
    std::ofstream fout;
    fout.open(name, std::ios::trunc);
    for (auto & i : data)
        fout << tuple_to_string(i) << '\n';
    fout.close();
}