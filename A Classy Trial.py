class Fruit(object):
	"""A class that makes various tasty fruits."""
	def __init__(self, name, color, flavor, poisonous):
		self.name = name
		self.colour = color
		self.flavor = flavor
		self.poisonous = poisonous
	
	def description(self):
		print ("I'm a {} {} and I taste {}.".format(self.colour, self.name, self.flavor))
		
	def is_edible(self):
		if not self.poisonous:
			print ("Yep! I'm edible.")
		else:
			print ("Don't eat me! I am super poisonous.")

lemon = Fruit("lemon", "yellow", "sour", False)

#lemon.description()
#lemon.is_edible()
print(lemon.colour)