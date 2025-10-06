function avgReward = evaluateAgent(agent, env, numEpisodes)
    rewards = zeros(numEpisodes,1);
    agent.AgentOptions.EpsilonGreedyExploration.Epsilon = 0;
    for i = 1:numEpisodes
        exp = sim(env, agent);
        rewards(i) = sum(exp.Reward);
    end
    avgReward = mean(rewards);
end
