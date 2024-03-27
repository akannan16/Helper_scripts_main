#Basic ggplot

ggplot(coverage_result_two_samples, aes(x=Sample_name, y=mean_coverage)) +
  geom_segment( aes(x=Sample_name, xend=Sample_name, y=0, yend=mean_coverage), color="grey", size = 1.5) +  labs(x = "Sample Name", y = "Mean Coverage") + geom_point( size=5, color="#69b3a2", fill=alpha("#69b3a2", 0.7), alpha=0.8, shape=20, stroke=2) +
  coord_flip() + theme_light()
theme(panel.grid.major.y = element_blank())
