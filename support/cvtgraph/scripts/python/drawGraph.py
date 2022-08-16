# This Python file uses the following encoding: utf-8
import plotly.graph_objects as go


# A graph drawer
class GraphDrawer:
    # Constructor
    def __init__(self, tit=''):
        self.title = tit
        self.node_x = []
        self.node_y = []
        self.node_color = []
        self.node_text = []
        self.node_hover_text = []
        self.edge_x = []
        self.edge_y = []

    # Add node
    def addNode(self, x, y, color=0, text='', hover_text=''):
        self.node_x.append(x)
        self.node_y.append(y)
        self.node_color.append(color)
        self.node_text.append(text)
        self.node_hover_text.append(hover_text)

    # add edge
    def addArc(self, x0, y0, x1, y1):
        self.edge_x.append(x0)
        self.edge_x.append(x1)
        self.edge_x.append(None)
        self.edge_y.append(y0)
        self.edge_y.append(y1)
        self.edge_y.append(None)

    # Draw
    def draw(self):
        edge_trace = go.Scatter(
            x=self.edge_x, y=self.edge_y,
            line=dict(width=0.5, color='#888'),
            hoverinfo='none',
            mode='lines')

        ##
        node_trace = go.Scatter(
            x=self.node_x, y=self.node_y,
            mode='text',
            hoverinfo='text'
            )

        ##
        node_trace.marker.color = self.node_color
        node_trace.text = self.node_text
        node_trace.hovertext = self.node_hover_text

        ##
        fig = go.Figure(data=[edge_trace, node_trace],
                        layout=go.Layout(
                        title=self.title,
                        titlefont_size=16,
                        showlegend=False,
                        hovermode='closest',
                        margin=dict(b=20,l=5,r=5,t=40),
                        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False))
                        )
        ##
        fig.update_layout(
            font_size=22,
            title_font_size=25,
            hoverlabel_font_size=32
        )
        ##
        fig.show()

