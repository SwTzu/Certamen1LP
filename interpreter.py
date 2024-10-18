import sys
from antlr4 import *
from parser.AsymmetricCALexer import AsymmetricCALexer
from parser.AsymmetricCAParser import AsymmetricCAParser
from parser.AsymmetricCAVisitor import AsymmetricCAVisitor
from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap

class AsymmetricCAInterpreter(AsymmetricCAVisitor):
    def __init__(self):
        self.cas = {}  # Stores CAs with their dimensions
        self.cells = defaultdict(dict)  # Stores cells within each CA
        self.connections = []  # Stores connections between cells
        self.parameters = {}  # Stores global parameters
        self.rules = ""  # Stores the transition rules
        self.time_unit = 1  # Default time unit
        self.simulation_history = []  # To store the state after each step

    def visitProgram(self, ctx):
        for stmt in ctx.statement():
            self.visit(stmt)
        return None

    def visitCaDeclaration(self, ctx):
        ca_name = ctx.ID().getText()
        dims = [int(dim.getText()) for dim in ctx.dimensions().INT()]
        self.cas[ca_name] = {'dimensions': dims}
        print(f"Defined CA '{ca_name}' with dimensions {dims}")
        return None

    def visitCellDeclaration(self, ctx):
        ca_name = ctx.ID().getText()
        pos = tuple(int(p.getText()) for p in ctx.cellPosition().INT())
        props = {}
        for prop in ctx.cellProperties().propertyAssignment():
            key = prop.ID().getText()
            value = prop.value().getText()
            props[key] = float(value) if '.' in value else int(value)
        self.cells[ca_name][pos] = props
        print(f"Defined cell in '{ca_name}' at {pos} with properties {props}")
        return None

    def visitConnection(self, ctx):
        ca1 = ctx.ID(0).getText()
        pos1 = tuple(int(p.getText()) for p in ctx.cellPosition(0).INT())
        ca2 = ctx.ID(1).getText()
        pos2 = tuple(int(p.getText()) for p in ctx.cellPosition(1).INT())
        self.connections.append(((ca1, pos1), (ca2, pos2)))
        print(f"Connected '{ca1}' at {pos1} with '{ca2}' at {pos2}")
        return None

    def visitParameterSetting(self, ctx):
        key = ctx.ID().getText()
        value = ctx.value().getText()
        self.parameters[key] = float(value) if '.' in value else int(value)
        print(f"Set parameter '{key}' to {self.parameters[key]}")
        return None

    def visitTimeUnitSetting(self, ctx):
        value = ctx.value().getText()
        self.time_unit = float(value) if '.' in value else int(value)
        print(f"Set time unit to {self.time_unit}")
        return None

    def visitAdvanceSimulation(self, ctx):
        if ctx.INT():
            steps = int(ctx.INT().getText())
        else:
            steps = 1  # Default to one step
        print(f"Advancing simulation by {steps} step(s).")
        self.advance_simulation(steps)
        return None

    def advance_simulation(self, steps):
        for step in range(steps):
            print(f"--- Simulation Step {step + 1} ---")
            self.simulate_step()
        # Visualize the simulation after all steps
        self.visualize_simulation()

    def simulate_step(self):
        # Create a copy of the cells to store updates
        new_cells = defaultdict(dict)

        # Apply transition rules and movement to each cell
        for ca_name, cells in self.cells.items():
            for pos, props in cells.items():
                # Initialize new cell properties if not already set
                if pos not in new_cells[ca_name]:
                    new_cells[ca_name][pos] = props.copy()

                # Retrieve compartment counts
                S = props.get('susceptible', 0)
                E = props.get('exposed', 0)
                I = props.get('infected', 0)
                R = props.get('recovered', 0)
                N = props.get('population', 1)  # Avoid division by zero

                # Disease parameters
                beta = self.parameters.get('infection_rate', 0.0)
                gamma = self.parameters.get('recovery_rate', 0.0)
                sigma = self.parameters.get('incubation_rate', 0.1)

                # Movement probability
                movement_prob = self.parameters.get('movement_probability', 0.0)

                # Compute new infections
                new_exposed = beta * S * I / N
                new_infected = sigma * E
                new_recovered = gamma * I

                # Update compartments
                S_new = S - new_exposed
                E_new = E + new_exposed - new_infected
                I_new = I + new_infected - new_recovered
                R_new = R + new_recovered

                # Ensure counts are non-negative
                S_new = max(S_new, 0)
                E_new = max(E_new, 0)
                I_new = max(I_new, 0)
                R_new = max(R_new, 0)

                # Movement between cells
                moved_S = moved_E = moved_I = moved_R = 0

                if movement_prob > 0:
                    # Get neighbors (now includes cross-CA connections)
                    neighbors = self.get_neighbors(ca_name, pos)

                    for neighbor_ca_name, neighbor_pos in neighbors:
                        neighbor_props = new_cells[neighbor_ca_name].get(neighbor_pos, self.cells[neighbor_ca_name].get(neighbor_pos, {}).copy())

                        # Move susceptible individuals
                        num_moving_S = int(S_new * movement_prob / len(neighbors))
                        if num_moving_S > 0:
                            S_new -= num_moving_S
                            neighbor_props['susceptible'] = neighbor_props.get('susceptible', 0) + num_moving_S
                            moved_S += num_moving_S

                        # Move exposed individuals
                        num_moving_E = int(E_new * movement_prob / len(neighbors))
                        if num_moving_E > 0:
                            E_new -= num_moving_E
                            neighbor_props['exposed'] = neighbor_props.get('exposed', 0) + num_moving_E
                            moved_E += num_moving_E

                        # Move infected individuals
                        num_moving_I = int(I_new * movement_prob / len(neighbors))
                        if num_moving_I > 0:
                            I_new -= num_moving_I
                            neighbor_props['infected'] = neighbor_props.get('infected', 0) + num_moving_I
                            moved_I += num_moving_I

                        # Move recovered individuals
                        num_moving_R = int(R_new * movement_prob / len(neighbors))
                        if num_moving_R > 0:
                            R_new -= num_moving_R
                            neighbor_props['recovered'] = neighbor_props.get('recovered', 0) + num_moving_R
                            moved_R += num_moving_R

                        # Update neighbor's properties
                        new_cells[neighbor_ca_name][neighbor_pos] = neighbor_props
                # Update current cell's properties
                new_cells[ca_name][pos]['susceptible'] = int(S_new)
                new_cells[ca_name][pos]['exposed'] = int(E_new)
                new_cells[ca_name][pos]['infected'] = int(I_new)
                new_cells[ca_name][pos]['recovered'] = int(R_new)

                # Print movement information
                print(f"Cell {ca_name}{pos} moved: S={moved_S}, E={moved_E}, I={moved_I}, R={moved_R}")
                print(f"Updated cell {ca_name}{pos}: {new_cells[ca_name][pos]}")

        # Update the cells with new values
        self.cells = new_cells

        # Save the current state after the step
        self.simulation_history.append(self.copy_cells_state())

    def copy_cells_state(self):
        # Create a deep copy of the cells to avoid reference issues
        state_copy = defaultdict(dict)
        for ca_name, cells in self.cells.items():
            for pos, props in cells.items():
                state_copy[ca_name][pos] = props.copy()
        return state_copy

    def get_neighbors(self, ca_name, pos):
        neighbors = []

        # All possible relative positions in a Moore neighborhood
        relative_positions = [(-1, -1), (-1, 0), (-1, 1),
                            (0, -1),          (0, 1),
                            (1, -1),  (1, 0),  (1, 1)]

        # Positions of all cells in the CA
        all_positions = set(self.cells[ca_name].keys())

        # Include connected cells from other CAs
        connected_cells = []
        for connection in self.connections:
            if (ca_name, pos) == connection[0]:
                connected_cells.append(connection[1])
            elif (ca_name, pos) == connection[1]:
                connected_cells.append(connection[0])

        # Check each relative position to find valid neighbors
        for dx, dy in relative_positions:
            neighbor_pos = (pos[0] + dx, pos[1] + dy)
            if neighbor_pos in all_positions:
                neighbors.append((ca_name, neighbor_pos))

        # Add connected cells as neighbors
        neighbors.extend(connected_cells)

        return neighbors

    def visualize_step(self, step_index):
        state = self.simulation_history[step_index]
        for ca_name in self.cas:
            # Get all positions defined in the CA
            positions = self.cells[ca_name].keys()
            xs, ys = zip(*positions)
            max_x, max_y = max(xs), max(ys)

            # Create a grid large enough to include all positions
            grid = np.full((max_x + 1, max_y + 1), np.nan)

            # Map disease states to integers for visualization
            state_mapping = {'susceptible': 0, 'exposed': 1, 'infected': 2, 'recovered': 3}

            for pos, props in state[ca_name].items():
                x, y = pos
                counts = {state: props.get(state, 0) for state in state_mapping.keys()}
                dominant_state = max(counts, key=counts.get)
                grid[x, y] = state_mapping[dominant_state]

            # Custom color map
            colors = ['blue', 'orange', 'red', 'green']
            cmap = ListedColormap(colors)

            plt.figure(figsize=(6, 6))
            plt.imshow(grid, cmap=cmap, origin='lower')

            # Overlay counts
            for (x, y), props in state[ca_name].items():
                counts_text = f"S:{int(props.get('susceptible', 0))}\nE:{int(props.get('exposed', 0))}\nI:{int(props.get('infected', 0))}\nR:{int(props.get('recovered', 0))}"
                plt.text(y, x, counts_text, ha='center', va='center', color='white', fontsize=8)

            plt.title(f"{ca_name} - Simulation Step {step_index + 1}")
            plt.colorbar(ticks=range(len(state_mapping)), label='Disease State', boundaries=np.arange(-0.5, len(state_mapping) + 0.5, 1))
            plt.clim(-0.5, len(state_mapping) - 0.5)
            plt.xticks(range(max_y + 1))
            plt.yticks(range(max_x + 1))
            plt.gca().set_xticks(np.arange(-0.5, max_y + 1, 1), minor=True)
            plt.gca().set_yticks(np.arange(-0.5, max_x + 1, 1), minor=True)
            plt.grid(which='minor', color='gray', linestyle='-', linewidth=1)
            plt.tight_layout()
            plt.show()


    def visualize_simulation(self):
        for step_index in range(len(self.simulation_history)):
            self.visualize_step(step_index)

    def visitRuleDefinition(self, ctx):
        # Extract the text between '{' and '}'
        rule_content = ctx.ruleBody().getText()
        self.rules = rule_content
        print("Defined transition rules.")
        return None

def main(argv):
    input_stream = FileStream(argv[1])
    lexer = AsymmetricCALexer(input_stream)
    stream = CommonTokenStream(lexer)
    parser = AsymmetricCAParser(stream)
    tree = parser.program()

    interpreter = AsymmetricCAInterpreter()
    interpreter.visit(tree)

    # For demonstration, print the internal state
    print("\nFinal State:")
    print("CAs:", interpreter.cas)
    print("Cells:", dict(interpreter.cells))
    print("Connections:", interpreter.connections)
    print("Parameters:", interpreter.parameters)
    print("Time Unit:", interpreter.time_unit)
    print("Rules:", interpreter.rules)

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python interpreter.py <filename>")
    else:
        main(sys.argv)
