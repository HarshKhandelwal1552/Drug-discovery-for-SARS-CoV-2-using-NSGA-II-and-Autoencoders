import random
from id import Individual
class Problem:

    def __init__(self, objectives, num_of_variables, expand=True, same_range=False):

        self.num_of_objectives = len(objectives)
        self.num_of_variables = num_of_variables
        self.objectives = objectives
        self.expand = expand


    def generate_individual(self, data):
        individual = Individual()
        individual.features = [data]
        return individual

    def calculate_objectives(self, individual):
        if self.expand:
            individual.objectives = [f(*individual.features) for f in self.objectives]
        else:
            individual.objectives = [f(individual.features) for f in self.objectives]