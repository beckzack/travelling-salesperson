#!/usr/bin/python3

from which_pyqt import PYQT_VER
if PYQT_VER == 'PYQT5':
	from PyQt5.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT4':
	from PyQt4.QtCore import QLineF, QPointF
else:
	raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))

import time
import numpy as np
from TSPClasses import *
import heapq
import itertools

class TSPSolver:
	def __init__( self, gui_view ):
		self._scenario = None

	def setupWithScenario( self, scenario ):
		self._scenario = scenario


	''' <summary>
		This is the entry point for the default solver
		which just finds a valid random tour.  Note this could be used to find your
		initial BSSF.
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of solution, 
		time spent to find solution, number of permutations tried during search, the 
		solution found, and three null values for fields not used for this 
		algorithm</returns> 
	'''
	
	def defaultRandomTour( self, time_allowance=60.0 ):
		results = {}
		cities = self._scenario.getCities()
		ncities = len(cities)
		foundTour = False
		count = 0
		bssf = None
		start_time = time.time()
		while not foundTour and time.time()-start_time < time_allowance:
			# create a random permutation
			perm = np.random.permutation( ncities )
			route = []
			# Now build the route using the random permutation
			for i in range( ncities ):
				route.append( cities[ perm[i] ] )
			bssf = TSPSolution(route)
			count += 1
			if bssf.cost < np.inf:
				# Found a valid route
				foundTour = True
		end_time = time.time()
		results['cost'] = bssf.cost if foundTour else math.inf
		results['time'] = end_time - start_time
		results['count'] = count
		results['soln'] = bssf
		results['max'] = None
		results['total'] = None
		results['pruned'] = None
		return results


	''' <summary>
		This is the entry point for the greedy solver, which you must implement for 
		the group project (but it is probably a good idea to just do it for the branch-and
		bound project as a way to get your feet wet).  Note this could be used to find your
		initial BSSF.
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of best solution, 
		time spent to find best solution, total number of solutions found, the best
		solution found, and three null values for fields not used for this 
		algorithm</returns> 
	'''

	def greedy( self,time_allowance=60.0 ):
		results = {}
		cities = self._scenario.getCities()
		ncities = len(cities)
		count = 0
		bssf = None
		start_time = time.time()

		while time.time() - start_time < time_allowance:
			for city in cities:
				route = []
				start_city = city
				# Use a set rather than a list for visited because sets have O(1) average search time complexity
				visited = set()
				route.append(city)
				visited.add(city)

				while len(route) < ncities:
					next_city = self.getNextCity(city, filter(lambda c: c not in visited, cities))
					if next_city == None:
						break
					route.append(next_city)
					visited.add(next_city)
					city = next_city

				# No route found
				if len(route) != ncities or city.costTo(start_city) == np.inf:
					continue

				temp_solution = TSPSolution(route)
				count += 1

				if bssf == None or bssf.cost > temp_solution.cost:
					bssf = temp_solution
				
			results['cost'] = bssf.cost
			results['time'] = time.time() - start_time
			results['count'] = count
			results['soln'] = bssf
			results['max'] = None
			results['total'] = None
			results['pruned'] = None
			
			return results
		
		results['cost'] = np.inf
		results['time'] = time.time() - start_time
		results['count'] = count
		results['soln'] = bssf
		results['max'] = None
		results['total'] = None
		results['pruned'] = None
		
		return results

	def getNextCity(self, currCity, cities):
		minCost = math.inf
		nextCity = None

		for city in cities:
			if city._name == currCity._name:
				continue
			
			if currCity.costTo(city) < minCost:
				nextCity = city
				minCost = currCity.costTo(city)
		
		return nextCity


	''' <summary>
		This is the entry point for the branch-and-bound algorithm that you will implement
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of best solution, 
		time spent to find best solution, total number solutions found during search (does
		not include the initial BSSF), the best solution found, and three more ints: 
		max queue size, total number of states created, and number of pruned states.</returns> 
	'''
		
	def branchAndBound( self, time_allowance=60.0 ):
		results = {}
		self.cities = self._scenario.getCities()
		ncities = len(self.cities)
		count = 0
		heap_queue = []
		max_queue_size = 0
		num_pruned = 0
		num_states = 0
		start_time = time.time()

		# Get the initial best solution so far using the greedy algorithm
		# The greedy algorithm is at worst O(n^2) time to find a solution where n is the number of cities
		bssf = self.greedy().get('soln')
		count += 1
		lower_bound = bssf.cost

		# Creating and populating the cost matrix takes O(n^2) time where n is the number of cities
		matrix = np.zeros((ncities, ncities))

		for i in range(ncities):
			for j in range(ncities):
				if i == j:
					matrix[i, j] = np.inf
					continue
				matrix[i, j] = self.cities[i].costTo(self.cities[j])
		
		initial_node = self.reduceMatrix(DataWrapperClass(matrix, 0, [], 0), 0, True)
		matrix = initial_node.matrix.copy()
		heap_queue = [initial_node]
		heapq.heapify(heap_queue)

		# Find the path using branch and bound
		# Worst case time complexity is O(n^2 * n!) but were hoping that our lower bound and bssf are good enough
		# to avoid that. Average time complexity is O(n^2 * b^n), where b is the average number of nodes put on the
		# queue and n is the overall number of nodes.
		while time.time() - start_time < time_allowance:
			if len(heap_queue) == 0:
				break
			if len(heap_queue) > max_queue_size:
				max_queue_size = len(heap_queue)

			# Now we pop the minimum from the heap, and from what I've found online it looks like they keep track
			# of the minimum index so the worst case time complexity is O(log(n))
			curr_node = heapq.heappop(heap_queue)
			for i in range(ncities):
				if i in curr_node.route:
					continue

				if curr_node.matrix[curr_node.index, i] != np.inf:
					num_states += 1
					comp_node = self.reduceMatrix(curr_node, i)

					if comp_node.cost < lower_bound:
						heapq.heappush(heap_queue, comp_node)

						# If we found a solution we know that it will be the bssf because the cost was less than the
						# lower bound
						if len(comp_node.route) == ncities:
							count += 1
							bssf = DataWrapperClass(comp_node.matrix.copy(), comp_node.cost, comp_node.route.copy(), comp_node.index)
					else:
						num_pruned += 1

		final_route = []
		for index in bssf.route:
			# I wish I could tell you why I need this check but I can't
			# The algorithm would work for some problems using the indeces but wouldn't work for other ones
			if type(index) == int:
				final_route.append(self.cities[index])
			else:
				final_route.append(index)
		
		end_time = time.time()
		results['cost'] = bssf.cost
		results['time'] = end_time - start_time
		results['count'] = count
		results['soln'] = TSPSolution(final_route)
		results['max'] = max_queue_size
		results['total'] = num_states
		results['pruned'] = num_pruned

		return results

	def reduceMatrix(self, data_object, dest_index, first=False):
		matrix = data_object.matrix.copy()
		start = data_object.index
		reduction_cost = data_object.cost
		
		#For the first matrix reduction we only subtract the minimum row and column values from the other elements
		# For every reduction after that we also set the rows and columns of our city pairing to infinity
		if not first:
			reduction_cost += matrix[start, dest_index]
			matrix[start] = np.inf
			matrix[:,dest_index] = np.inf
			matrix[dest_index, start] = np.inf

		# Find the minimum in each column and row
		for i in range(len(matrix)):
			row_min = matrix[i].min()
			if row_min != np.inf:
				reduction_cost += row_min
				matrix[i] = matrix[i] - row_min

		for i in range(len(matrix)):
			col_min = matrix[:,i].min()
			if col_min != np.inf:
				reduction_cost += col_min
				matrix[:,i] = matrix[:,i] - col_min

		route = data_object.route.copy()
		route.append(dest_index)

		return DataWrapperClass(matrix, reduction_cost, route, dest_index)

	''' <summary>
		This is the entry point for the algorithm you'll write for your group project.
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of best solution, 
		time spent to find best solution, total number of solutions found during search, the 
		best solution found.  You may use the other three field however you like.
		algorithm</returns> 
	'''
		
	def fancy( self,time_allowance=60.0 ):
		pass
		
class DataWrapperClass():
	def __init__(self, matrix, cost, route, index):
		self.matrix = matrix
		self.cost = cost
		self.route = route
		self.index = index

	# Define a less than function so the heap knows how to organize itself
	def __lt__(self, comp):
		# We want to find all states on one level before moving to the next one
		# So we compare using the length of the partial paths as well as the cost
		if len(self.route) > len(comp.route):
			return True
		if len(self.route) == len(comp.route) and self.cost < comp.cost:
			return True
		return False
