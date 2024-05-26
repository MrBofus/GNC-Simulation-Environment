import getopt, sys

argumentList = sys.argv[1:]

print(argumentList)

# parser = argparse.ArgumentParser()

# parser.add_argument('--run-visualizer', dest='run_visualizer', action='store_true', default=False)
# parser.add_argument('--orbit', dest='orbit', action='store_true', default=False)

# args = parser.parse_args()

if not '--orbit' in argumentList:

    import software_example_simulation
    
    if '--run-visualizer' in argumentList:
        _ = input("\n\nRun visualizer? ")
        import visualizer

else:

    import orbit_modular_example

    for arg in argumentList:
        if '=' in arg:
            kwarglist.append(arg.split)
    # orbit_modular_example.main()

    if '--run-visualizer' in argumentList:
        _ = input("\n\nRun visualizer? ")
        import visualizer