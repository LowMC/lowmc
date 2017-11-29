#!/usr/bin/python3
''' This script determines the number of rounds needed for security for
    a given set of LowMC parameters. Use with '-h' flag for more information.
'''

from math import log2, ceil
import argparse
import itertools


# Define what to consider a negligible probability
NEGL_PROB = 2**(-100)


def main():
    ''' Determine the rounds needed for security for the given
        parameters.
    '''
    params = parse_program_arguments()
    check_parameter_validity(params)
    if params.verbosity == 'verbose':
        print('-' * 46)
        print('LowMC rounds determination')
        print('-' * 46)
        params.print()
        print('-' * 46)

    if params.verbosity == 'verbose':
        print('Calculating statistical rounds')
    statistical_rounds = determine_statistical_rounds(params)
    if params.verbosity == 'verbose':
        print('Calculating boomerang rounds')
    boomerang_rounds = determine_boomerang_rounds(params)
    if params.verbosity == 'verbose':
        print('Calculating derivative rounds')
    derivative_rounds = determine_derivative_rounds(params)
    if params.verbosity == 'verbose':
        print('Calculating interpolation rounds')
    interpolation_rounds = determine_interpolation_rounds(params)
    if params.verbosity == 'verbose':
        print('Calculating round-key guessing rounds')
    keyguess_state_rounds = determine_keyguess_state_rounds(params)
    keyguess_bit_rounds = determine_keyguess_bit_rounds(params, 1)
    if params.verbosity == 'verbose':
        print('Calculating polytopic attack rounds')
    polytopic_rounds = determine_polytopic_attack_rounds(params)

    distinguishers = []
    distinguishers.append(('Statistical with state guessing',
                           statistical_rounds + keyguess_state_rounds))
    distinguishers.append(('Boomerang attack',
                           boomerang_rounds))
    distinguishers.append(('Derivative + bit guessing',
                           derivative_rounds + keyguess_bit_rounds))
    distinguishers.append(('Derivative + interpolation',
                           derivative_rounds + interpolation_rounds))
    distinguishers.append(('Impossible polytopic attack',
                           polytopic_rounds))

    print_rounds(params, distinguishers)


###########################################################
# Determining secure rounds for a range of distinguishers
###########################################################
def determine_statistical_rounds(params):
    ''' Determine the number of rounds for which no good differential
        or linear trail exists with high probability.
    '''
    # Determine the number of rounds using a divide-and-conquer approach
    for round_exponent in itertools.count(0):
        if no_good_trail_after_round(params, 2**round_exponent):
            upper_bound = 2**round_exponent
            lower_excl_bound = 2**(round_exponent - 1)
            break
    while lower_excl_bound + 1 < upper_bound:
        rounds = lower_excl_bound + (upper_bound - lower_excl_bound) // 2
        if no_good_trail_after_round(params, rounds):
            upper_bound = rounds
        else:
            lower_excl_bound = rounds
    return upper_bound


def determine_boomerang_rounds(params):
    ''' Determine the number of rounds for which no good boomerang
        consisting of one differential trail for the top part and one
        for the bottom part can be constructed.
    '''
    # Determine the number of rounds using a divide-and-conquer approach
    for round_exponent in itertools.count(1):
        if no_good_boomerang_after_round(params, 2**round_exponent):
            upper_bound = 2**round_exponent
            lower_excl_bound = 2**(round_exponent - 1)
            break
    while lower_excl_bound + 1 < upper_bound:
        rounds = lower_excl_bound + (upper_bound - lower_excl_bound) // 2
        if no_good_boomerang_after_round(params, rounds):
            upper_bound = rounds
        else:
            lower_excl_bound = rounds
    return upper_bound


def determine_impossible_rounds(params):
    ''' Determine the number of rounds after which it should not be possible
        to find an impossible differential.
    '''
    rounds = determine_free_rounds(params, 1)
    rounds += 2 * polytopic_listing_rounds(params, 1)
    rounds += min(determine_free_rounds(params, 1),\
                  ceil((params.data_complexity - 1) / (3 * params.sboxes)))
    return rounds


def determine_derivative_rounds(params):
    ''' Determines the number of rounds after which there should be no
        derivative that always evaluates to zero.
    '''
    degree_rounds = determine_degree_rounds(params)
    influence_rounds = determine_influence_rounds(params)
    return degree_rounds + influence_rounds


def determine_polytopic_attack_rounds(params):
    ''' Determine the number of rounds after which an impossible
        polytopic attack is impossible.
    '''
    attacked_rounds = []
    for ddiff_size in range(1, ceil(2 * params.keysize / params.blocksize) + 1):
        if log2(ddiff_size + 1) > params.data_complexity:
            continue # Initial free rounds
        rounds = determine_free_rounds(params, ceil(log2(ddiff_size + 1)))
        # d-difference diffusion rounds
        rounds += polytopic_listing_rounds(params, ddiff_size)
        # Backwards key-guessing rounds
        rounds += polytopic_listing_rounds(params, ddiff_size) \
                   + determine_free_rounds(params, params.blocksize \
                   - params.data_complexity // ddiff_size)
        attacked_rounds.append(rounds)
    return max(attacked_rounds)


def determine_interpolation_rounds(params):
    ''' Determine the number of rounds needed for security against
        interpolation attacks. This is done by testing whether
        solving the linear system of equations takes more time
        or data than allowed.
    '''
    for rounds in itertools.count(1):
        terms = interpolation_terms(params, rounds)
        if log2(terms) >= params.keysize/2.3 \
           or log2(terms) >= params.data_complexity:
        # the 2.3 is used here as a lower bound on the complexity of solving
        # a system of linear equations.
            return rounds


def determine_keyguess_state_rounds(params):
    ''' Determine how many rounds back it is possible to guess the full
        state given the computationals constraints.
    '''
    return params.keysize // (3 * params.sboxes)


def determine_keyguess_bit_rounds(params, dimension):
    ''' Determine how many rounds back it is possible to guess a subspace
        of dimension 'dimension' of the state given the computationals
        constraints.
    '''
    free_rounds = determine_free_rounds(params, dimension)
    guess_state_rounds = determine_keyguess_state_rounds(params)
    return free_rounds + guess_state_rounds


####################################################
# Statistical and boomerang distinguisher functions
####################################################
def no_good_trail_after_round(params, rounds):
    ''' Determines whether a good trail exists after 'rounds' rounds.
    '''
    # A "good" trail should have a probability higher than
    # 2^-(data_complexity)
    max_active_sboxes = params.data_complexity // 2
    all_good_trails = all_possible_good_trails(params, max_active_sboxes, rounds)
    # To realize a trail the linear layers have to connect corresponding
    # differences. The inverse of the probability that this happens for
    # any given trail consisting of 'rounds' single round differentials
    # is now calculated as:
    inv_realization_probability = (2**params.blocksize - 1)**(rounds - 1)
    # We would like that the probability that at least one trail is
    # realized times the number of all good trails is smaller than the
    # negligible probability.
    return int(1/NEGL_PROB) * all_good_trails < inv_realization_probability


def all_possible_good_trails(params, max_active_sboxes, rounds):
    ''' Returns an upper bound for the probability that
        there exist a trails over 'rounds' rounds that
        activates at most 'max_active_sboxes' Sboxes.
    '''
    current_trails = [0 for _ in range(max_active_sboxes + 1)]
    # Store the number of good trails after 1 round
    for active_sboxes in range(max_active_sboxes + 1):
        current_trails[active_sboxes] = one_round_trails(params, active_sboxes)
    for _ in range(2, rounds + 1):
        new_trails = [0 for _ in range(max_active_sboxes + 1)]
        for prev_actives in range(max_active_sboxes + 1):
            for new_actives in range(max_active_sboxes - prev_actives + 1):
                new_trails[prev_actives + new_actives] \
                  += current_trails[prev_actives] \
                     * one_round_trails(params, new_actives)
        current_trails = new_trails
    # Count the total number of valid trails
    all_good_trails = sum(current_trails)
    return all_good_trails


def no_good_boomerang_after_round(params, rounds):
    ''' Determines whether a good boomerang trail exists after 'rounds'
        rounds.
    '''
    # Since boomerang trails are used twice, each S-box contributes now
    # a probability of 2^-4.
    max_actives = params.data_complexity // 4
    top_rounds = rounds // 2
    bottom_rounds = rounds - top_rounds
    # Note: Boomerang probability for a sub-characteristic
    #       has to be a fourth of the normal bound
    combination_is_secure = []
    for top_actives in range(max_actives + 1):
        bottom_actives = max_actives - top_actives
        top_good_trails = all_possible_good_trails(params, top_actives, top_rounds)
        bottom_good_trails = \
                all_possible_good_trails(params, bottom_actives, bottom_rounds)
        inv_top_realization_prob = (2**params.blocksize - 1)**(top_rounds - 1)
        inv_bottom_realization_prob = (2**params.blocksize - 1)**(bottom_rounds - 1)
        if int(1/NEGL_PROB) * top_good_trails < inv_top_realization_prob or \
           int(1/NEGL_PROB) * bottom_good_trails < inv_bottom_realization_prob:
            combination_is_secure.append(True)
        else:
            combination_is_secure.append(False)
            break
    return all(combination_is_secure)


def one_round_trails(params, active_sboxes):
    ''' Number of one-round trails activating exactly 'active_sboxes'
        S-boxes.
    '''
    # There are exactly 4 possible output differences for a fixed input
    # to an S-box
    return activating_vectors(params, active_sboxes) * 4 ** active_sboxes


def activating_vectors(params, active_sboxes):
    ''' Number of vectors that activate the given number of active
        S-boxes out of all S-boxes with the given length of
        the identity part of the nonlinear layer
    '''
    return choose(params.sboxes, active_sboxes) * 7**active_sboxes \
           * 2**params.identity_bits


####################################################
# Impossible distinguisher functions
####################################################
def determine_free_rounds(params, dimension):
    ''' Determine the maximal number of rounds for which there exists a
        subspace of dimension 'dimension' that depends only linearly on
        the input bits.
    '''
    if dimension > params.blocksize:
        raise ValueError('dimension must not be larger than blocksize')
    if 3 * params.sboxes == params.blocksize:
        return 0
    else:
        # Add 1 as security margin
        return (params.blocksize - dimension) // (3 * params.sboxes) + 1


####################################################
# Polytopic distinguisher functions
####################################################
def polytopic_listing_rounds(params, ddiff_size):
    ''' Determine the number of rounds after which all d-differences are
        reachable or after which it should not be possible
        to list all reachable d-differences of size 'ddiff_size', given
        a fixed, random input d-difference.
    '''
    rounds = 0
    diffusion_per_round = \
            calculate_average_polytopic_diffusion(params, ddiff_size)
    diffusion = 0.0
    while diffusion < params.keysize and \
          diffusion < ddiff_size * params.blocksize:
        rounds += 1
        diffusion += diffusion_per_round
    return rounds

def calculate_average_polytopic_diffusion(params, ddiff_size):
    ''' Determine the average diffusion of d-difference of size
        'ddiff_size'.
    '''
    sboxes = params.sboxes
    # In the following, we calculate the number of new d-difference
    # generated by all possible patterns and divide by the number
    # of all possible patterns
    all_created_differences = 0
    # Number of S-box d-difference with at least two non-zero differences
    sbox_two_active = 8**ddiff_size - 1 - ddiff_size * 7
    for inactive in range(sboxes + 1):
        for one_active in range(sboxes - inactive + 1):
            number_of_patterns = choose(sboxes, inactive) * \
                    choose(sboxes - inactive, one_active) * \
                    (ddiff_size * 7)**one_active * \
                    sbox_two_active**(sboxes - inactive - one_active)
            new_differences = number_of_patterns * \
                    4**one_active * 8**(sboxes - inactive - one_active)
            all_created_differences += new_differences
    # Divide by number of all patterns and take log2
    all_created_differences = log2(all_created_differences)
    all_created_differences -= 3 * sboxes * ddiff_size
    return all_created_differences


####################################################
# Derivative distinguisher functions
####################################################
def determine_degree_rounds(params):
    ''' Determine the number of rounds needed so that the maximal possible
        degree is not smaller than the allowed data complexity minus one.
    '''
    for rounds in itertools.count(1):
        max_degree = determine_degree_upper_bound(params, rounds)
        if max_degree >= params.data_complexity - 1:
            return rounds


def determine_degree_upper_bound(params, rounds):
    ''' Calculates an upper bound for the algebraic degree
        after r rounds. m is the number of Sboxes per Sbox
        layer, n ist the block size.
    '''
    degree = 1
    for _ in range(rounds):
        degree = min(2*degree, params.sboxes+degree, \
                     (params.blocksize+degree)//2)
    return degree


def determine_influence_rounds(params):
    ''' Estimate number of rounds for one bit to influence all others
    '''
    return ceil(params.blocksize / (7 / 8 * params.sboxes * 3))


####################################################
# Interpolation distinguisher functions
####################################################
def interpolation_terms(params, rounds):
    ''' Estimate the number of different terms in the key bits that appear
        in an interpolation attack on the last 'rounds' rounds
    '''
    keybit_terms = [0 for _ in range(params.blocksize + 1)]
    # keybit_terms[d] will hold an estimation of the number of terms
    # in the key bits of degree d
    # After 1 round
    keybit_terms[0] = 1
    keybit_terms[1] = params.blocksize
    keybit_terms[2] = 3 * params.sboxes
    # For each additional round
    for _ in range(1, rounds):
        # Store the new number of terms in newkeybit_terms
        newkeybit_terms = []
        # The number of constant terms and linear terms is always 1 and n
        newkeybit_terms.append(1)
        newkeybit_terms.append(params.blocksize)
        # The terms of higher degree are generated by multiplying terms
        # of lower degree
        for degree in range(2, params.blocksize + 1):
            terms_of_degree = 0
            for degree_1st_factor in range(0, degree // 2 + 1):
                terms_of_degree += keybit_terms[degree_1st_factor] \
                                   * keybit_terms[degree - degree_1st_factor]
            newkeybit_terms.append(
                min(terms_of_degree, choose(params.blocksize, degree)))
            # The number of terms of degree 'degree' is always upper bounded
            # by choose(blocksize, degree).
        keybit_terms = newkeybit_terms
    # To estimate the number of terms we combine the estimate
    # for the number of terms in the ciphertext bits with the number
    # of terms possible for the key bits.
    # As a term of degree 'degree' in the ciphertext bits can only have a
    # coefficient of degree 2**rounds-degree in the key bits, we have the
    # following
    terms = 0
    for degree in range(min(2**rounds + 1, params.blocksize + 1)):
        terms += min(keybit_terms[degree], terms_with_bounded_degree(params.keysize, \
                   2**rounds - degree))
    return terms


def terms_with_bounded_degree(variables, max_degree):
    ''' Calculate the number of possible terms of with the given maximal
        degree in the given number of variables.
    '''
    terms = 0
    for degree in range(max_degree + 1):
        terms += choose(variables, degree)
    return terms


###################################
# Helper functions
###################################
def parse_program_arguments():
    ''' Take the command line arguments and attribute them to the
        corresponding variables.
    '''
    program_description = 'Calculate the number of rounds for LowMC' \
                          'in dependence on a given parameter set.'
    parser = argparse.ArgumentParser(description=program_description)
    parser.add_argument('block', type=int,
                        help='specifies block size in bits')
    parser.add_argument('sboxes', type=int,
                        help='specifies number of Sboxes per nonlinear layer')
    parser.add_argument('data', type=int,
                        help='specifies the log2 of allowed data complexity')
    parser.add_argument('key', type=int,
                        help='specifies key size in bits')
    verbosity_group = parser.add_mutually_exclusive_group()
    verbosity_group.add_argument('-v', '--verbose', action='store_true',
                                 help='print additional information')
    verbosity_group.add_argument('-q', '--quiet', action='store_true',
                                 help='only print the number of total rounds')
    parser.add_argument('--multiplicity', '-m',\
                        help='print multiplicative complexities',\
                        action='store_true')
    arguments = parser.parse_args()
    params = Parameters(blocksize=arguments.block,
                        sboxes=arguments.sboxes,
                        data_complexity=arguments.data,
                        keysize=arguments.key)
    if arguments.verbose:
        params.verbosity = 'verbose'
    elif arguments.quiet:
        params.verbosity = 'quiet'
    else:
        params.verbosity = 'normal'
    if arguments.multiplicity:
        params.with_multiplicity = True
    return params


def check_parameter_validity(params):
    ''' Check if the given parameters are coherent and valid.
    '''
    if params.blocksize < params.data_complexity \
       or params.sboxes * 3 > params.blocksize \
       or params.data_complexity > params.keysize \
       or any(p < 1 for p in {params.sboxes, params.blocksize,\
                              params.data_complexity, params.keysize}):
        print("Invalid parameter set")
        exit()


def print_rounds(params, distinguishers):
    ''' Print information the total rounds for the parameter set, the
        rounds for each distinguisher and the parameter set, dependent
        on the verbosity level.
    '''
    total_rounds = max(d[1] for d in distinguishers)
    if params.verbosity != 'quiet':
        print('-' * 46)
        print('Distinguisher', ' ' * 27, 'Rounds', sep='')
        print('-' * 46)
        for dist in distinguishers:
            print(dist[0], ' ' * (40 - len(dist[0])), sep='', end='')
            print('{:>6}'.format(dist[1]))
        print('-' * 46)
        print('Secure rounds:', ' ' * 26, sep='', end='')
        print('{:>6}'.format(total_rounds))
    else:
        print(total_rounds)
    if params.with_multiplicity:
        print('-' * 46)
        print('Total number of ANDs:', ' ' * 19, sep='', end='')
        print('{:>6}'.format(total_rounds * 3 * params.sboxes))
        print('Number of ANDs per bit:', ' ' * 17, sep='', end='')
        print('{:6.2f}'.format(total_rounds * 3.0 * params.sboxes \
                               / params.blocksize))
        print('AND-depth:', ' ' * 30, sep='', end='')
        print('{:>6}'.format(total_rounds))


def choose(z, n):
    ''' Calculates the binomial coefficient  "z choose n".
    '''
    numerator, denominator = 1, 1
    for i in range(n):
        numerator *= z - i
        denominator *= i + 1
    return numerator // denominator

#################################
# Memoization
#################################

class Memorizer:
    def __init__(self, function):
        self.f = function
        self.store = {}
    def __call__(self, *params):
        if params not in self.store:
            self.store[params] = self.f(*params)
        return self.store[params]

choose = Memorizer(choose)
one_round_trails = Memorizer(one_round_trails)
all_possible_good_trails = Memorizer(all_possible_good_trails)

#################################
# Parameter class
#################################

class Parameters:
    ''' This class contains the parameter set of a LowMC instantiation
    '''
    def __init__(self, blocksize, sboxes, data_complexity, keysize):
        self.blocksize = blocksize
        self.sboxes = sboxes
        self.data_complexity = data_complexity
        self.keysize = keysize
        self.identity_bits = blocksize - 3 * sboxes
        self.verbosity = 'normal'
        self.with_multiplicity = False
    def print(self):
        ''' Print the parameters.
        '''
        print('Block size:     ', self.blocksize)
        print('# of Sboxes:    ', self.sboxes)
        print('Data complexity:', self.data_complexity)
        print('Key size:       ', self.keysize)



##################################
# Executed code
##################################
if __name__ == '__main__':
    main()
