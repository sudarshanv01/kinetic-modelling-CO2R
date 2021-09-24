from catmap import ReactionModel 
mkm_file = 'aiida.mkm' 
model = ReactionModel(setup_file=mkm_file) 
model.output_variables += ['production_rate'] 
model.run() 
