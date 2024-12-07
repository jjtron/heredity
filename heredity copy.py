import csv
import itertools
import sys
import pprint

PROBS = {

    # Unconditional probabilities for having gene
    "gene": {
        2: 0.01,
        1: 0.03,
        0: 0.96
    },

    "trait": {

        # Probability of trait given two copies of gene
        2: {
            True: 0.65,
            False: 0.35
        },

        # Probability of trait given one copy of gene
        1: {
            True: 0.56,
            False: 0.44
        },

        # Probability of trait given no gene
        0: {
            True: 0.01,
            False: 0.99
        }
    },

    # Mutation probability
    "mutation": 0.01
}


def main():

    # Check for proper usage
    if len(sys.argv) != 2:
        sys.exit("Usage: python heredity.py data.csv")
    people = load_data(sys.argv[1])
    people = {
                'Harry': {'name': 'Harry', 'mother': 'Lily', 'father': 'James', 'trait': None},
                'James': {'name': 'James', 'mother': None, 'father': None, 'trait': None},
                'Lily': {'name': 'Lily', 'mother': None, 'father': None, 'trait': None}
            }
    '''
    SHOULD GET THESE OUTPUTS
        Harry 1  *  0.9503999999999999  *  0.9801
        Harry 0.9314870399999999
        James 0.9314870399999999  *  0.9503999999999999  *  0.96
        James 0.8498738715033598
        Lily 0.8498738715033598  *  0.9503999999999999  *  0.96
        Lily 0.7754113223777214
        
        :(  joint_probability returns correct results for no gene or trait in simple family
            expected joint probability to be in range [0.8664, 0.8864], got 0.7754113223777214 instead
    '''
    # Keep track of gene and trait probabilities for each person
    probabilities = {
        person: {
            "gene": {
                2: 0,
                1: 0,
                0: 0
            },
            "trait": {
                True: 0,
                False: 0
            }
        }
        for person in people
    }

    # Loop over all sets of people who might have the trait
    names = set(people)
    for have_trait in powerset(names):

        # Check if current set of people violates known information
        fails_evidence = any(
            (people[person]["trait"] is not None and
             people[person]["trait"] != (person in have_trait))
            for person in names
        )
        if fails_evidence:
            continue

        # Loop over all sets of people who might have the gene
        for one_gene in powerset(names):
            for two_genes in powerset(names - one_gene):

                # Update probabilities with new joint probability
                p = joint_probability(people, one_gene, two_genes, have_trait)
                update(probabilities, one_gene, two_genes, have_trait, p)

    # Ensure probabilities sum to 1
    normalize(probabilities)

    # Print results
    for person in people:
        print(f"{person}:")
        for field in probabilities[person]:
            print(f"  {field.capitalize()}:")
            for value in probabilities[person][field]:
                p = probabilities[person][field][value]
                print(f"    {value}: {p:.4f}")


def load_data(filename):
    """
    Load gene and trait data from a file into a dictionary.
    File assumed to be a CSV containing fields name, mother, father, trait.
    mother, father must both be blank, or both be valid names in the CSV.
    trait should be 0 or 1 if trait is known, blank otherwise.
    """
    data = dict()
    with open(filename) as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = row["name"]
            data[name] = {
                "name": name,
                "mother": row["mother"] or None,
                "father": row["father"] or None,
                "trait": (True if row["trait"] == "1" else
                          False if row["trait"] == "0" else None)
            }
    return data


def powerset(s):
    """
    Return a list of all possible subsets of set s.
    """
    s = list(s)
    return [
        set(s) for s in itertools.chain.from_iterable(
            itertools.combinations(s, r) for r in range(len(s) + 1)
        )
    ]


def joint_probability(people, one_gene, two_genes, have_trait):
    '''
    pprint.pp(people)
    pprint.pp(one_gene)
    pprint.pp(two_genes)
    pprint.pp(have_trait)
    '''
    """
    Compute and return a joint probability.

    The probability returned should be the probability that
        * everyone in set `one_gene` has one copy of the gene, and
        * everyone in set `two_genes` has two copies of the gene, and
        * everyone not in `one_gene` or `two_gene` does not have the gene, and
        * everyone in set `have_trait` has the trait, and
        * everyone not in set` have_trait` does not have the trait.
    """
    running_total = 1
    for person in people:
        fg = 0
        ft = 0
        if person not in one_gene and person not in two_genes:
            #print("*0",person)
            fg = calc_prob_gene(0, person in have_trait)
            ft = mutation_probabilities_zero_genes(people, person, one_gene, two_genes)
        if person in one_gene:
            #print("*1",person)
            fg = calc_prob_gene(1, person in have_trait)
            ft = mutation_probabilities_one_gene(people, person, one_gene, two_genes)
        if person in two_genes:
            #print("*2",person)
            fg = calc_prob_gene(2, person in have_trait)
            ft = mutation_probabilities_two_genes(people, person, one_gene, two_genes)
        #print(person, running_total, " * ", fg, " * ", ft)
        running_total = running_total * fg * ft
        #print(person, running_total)
    return running_total

def calc_prob_gene(gene_count, has_trait):
        #print(PROBS["gene"][gene_count], PROBS["trait"][gene_count][has_trait])
        return PROBS["gene"][gene_count] * PROBS["trait"][gene_count][has_trait]
        
def mutation_probabilities_zero_genes(people, the_person, one_gene, two_genes):
    mutation_probability = 0

    # Check if person has parents
    parent_list = parents(people, the_person)
    mother = parent_list["mother"]
    father = parent_list["father"]

    #print(the_person, "ZERO GENES ---> mother: ", mother, "; father: ", father)

    if mother == None and father == None:
        return PROBS['gene'][0]
        #print("No parents")

    if mother != None and father == None:
        print()
        #print("One parent")

    # Case: person has two parents
    if mother != None and father != None:
        #print("mother != None and father != None")
        # Case 1: one parent 0,0, other parent 1,1
        if mother not in one_gene and mother not in two_genes:
            if father in two_genes:
                mutation_probability = .0099
        if father not in one_gene and father not in two_genes:
            if mother in two_genes:
                mutation_probability = .0099

        # Case 2: mother 0,1 and father 0,1
        if mother in one_gene and father in one_gene:
            mutation_probability = .9901

        # Case 3: mother 0,0 and father 0,0
        if mother not in one_gene and mother not in two_genes:
            if father not in one_gene and father not in two_genes:
                mutation_probability = .9801

        # Case 4: mother 1,1 and father 1,1
        if mother in two_genes and father in two_genes:
            mutation_probability = .0001

        # Case 5: one parent 0,1, other parent 1,1
        if (mother in one_gene and father in two_genes) or (father in one_gene and mother in two_genes):
            mutation_probability = .01

        # Case 6: one parent 0,1, other parent 0,0
        if mother not in one_gene and mother not in two_genes:
            if father in one_gene:
                mutation_probability = .9900
        if father not in one_gene and father not in two_genes:
            if mother in one_gene:
                mutation_probability = .9900

    return mutation_probability * PROBS["trait"][0][False]


def mutation_probabilities_one_gene(people, the_person, one_gene, two_genes):
    mutation_probability = 0

    # Check if person has parents
    parent_list = parents(people, the_person)
    mother = parent_list["mother"]
    father = parent_list["father"]

    #print(the_person, "ONE GENE ---> mother: ", mother, "; father: ", father)

    if mother == None and father == None:
        return PROBS['gene'][1]
        #print("No parents")

    if mother != None and father == None:
        print()
        #print("One parent")

    # Case: person has two parents
    if mother != None and father != None:
        #print("mother != None and father != None")
        # Case 1: one parent 0,0, other parent 1,1
        if mother not in one_gene and mother not in two_genes:
            if father in two_genes:
                mutation_probability = .9802
        if father not in one_gene and father not in two_genes:
            if mother in two_genes:
                mutation_probability = .9802

        # Case 2: mother 0,1 and father 0,1
        if mother in one_gene and father in one_gene:
            mutation_probability = .9999

        # Case 3: mother 0,0 and father 0,0
        if mother not in one_gene and mother not in two_genes:
            if father not in one_gene and father not in two_genes:
                mutation_probability = .0099

        # Case 4: mother 1,1 and father 1,1
        if mother in two_genes and father in two_genes:
            mutation_probability = .0099

        # Case 5: one parent 0,1, other parent 1,1
        if (mother in one_gene and father in two_genes) or (father in one_gene and mother in two_genes):
            mutation_probability = .01

        # Case 6: one parent 0,1, other parent 0,0
        if mother not in one_gene and mother not in two_genes:
            if father in one_gene:
                mutation_probability = .9900
        if father not in one_gene and father not in two_genes:
            if mother in one_gene:
                mutation_probability = .9900

    return mutation_probability * PROBS["trait"][1][False]

def mutation_probabilities_two_genes(people, the_person, one_gene, two_genes):
    mutation_probability = 0

    # Check if person has parents
    parent_list = parents(people, the_person)
    mother = parent_list["mother"]
    father = parent_list["father"]

    #print(the_person, "TWO GENES ---> mother: ", mother, "; father: ", father)

    if mother == None and father == None:
        return PROBS['gene'][2]
        #print("No parents")

    if mother != None and father == None:
        print()
        #print("One parent")

    # Case: person has two parents
    if mother != None and father != None:
        #print("mother != None and father != None")
        # Case 1: one parent 0,0, other parent 1,1
        if mother not in one_gene and mother not in two_genes:
            if father in two_genes:
                mutation_probability = .0099
        if father not in one_gene and father not in two_genes:
            if mother in two_genes:
                mutation_probability = .0099

        # Case 2: mother 0,1 and father 0,1
        if mother in one_gene and father in one_gene:
            mutation_probability = .991

        # Case 3: mother 0,0 and father 0,0
        if mother not in one_gene and mother not in two_genes:
            if father not in one_gene and father not in two_genes:
                mutation_probability = .001

        # Case 4: mother 1,1 and father 1,1
        if mother in two_genes and father in two_genes:
            mutation_probability = .9801

        # Case 5: one parent 0,1, other parent 1,1
        if (mother in one_gene and father in two_genes) or (father in one_gene and mother in two_genes):
            mutation_probability = .99

        # Case 6: one parent 0,1, other parent 0,0
        if mother not in one_gene and mother not in two_genes:
            if father in one_gene:
                mutation_probability = .01
        if father not in one_gene and father not in two_genes:
            if mother in one_gene:
                mutation_probability = .01

    return mutation_probability * PROBS["trait"][2][False]


def parents(people, person):
    # Check if person has parents
    retval = {}
    retval["mother"] = people[person]["mother"]
    retval["father"] = people[person]["father"]
    return retval


def update(probabilities, one_gene, two_genes, have_trait, p):
    '''
    pprint.pp(probabilities)
    pprint.pp(one_gene)
    pprint.pp(two_genes)
    pprint.pp(have_trait)
    pprint.pp(p)
    '''
    """
    Add to `probabilities` a new joint probability `p`.
    Each person should have their "gene" and "trait" distributions updated.
    Which value for each distribution is updated depends on whether
    the person is in `have_gene` and `have_trait`, respectively.
    """
    for person in probabilities:
        if person not in one_gene and person not in two_genes:
            probabilities[person]['gene'][0] += p
        if person in one_gene:
            probabilities[person]['gene'][1] += p
        if person in two_genes:
            probabilities[person]['gene'][2] += p

        if person in have_trait:
            probabilities[person]['trait'][True] += p
        else:
            probabilities[person]['trait'][False] += p


def normalize(probabilities):

    """
    Update `probabilities` such that each probability distribution
    is normalized (i.e., sums to 1, with relative proportions the same).
    """
    for person in probabilities:
        alpha = 1 / (probabilities[person]['gene'][0] + probabilities[person]['gene'][1] + probabilities[person]['gene'][2])
        probabilities[person]['gene'][0] = alpha * probabilities[person]['gene'][0]
        probabilities[person]['gene'][1] = alpha * probabilities[person]['gene'][1]
        probabilities[person]['gene'][2] = alpha * probabilities[person]['gene'][2]

        alpha = 1 / (probabilities[person]['trait'][True] + probabilities[person]['trait'][False])
        probabilities[person]['trait'][True] = alpha * probabilities[person]['trait'][True]
        probabilities[person]['trait'][False] = alpha * probabilities[person]['trait'][False]


if __name__ == "__main__":
    main()
