import os
from mkdocs.plugins import BasePlugin
# from mkdocs.structure.pages import MarkdownPage


class ArgparseToMDPlugin(BasePlugin):

    def on_page_read_source(self, page, config):
        # if not isinstance(page, MarkdownPage):
        # return None
        content = None
        # Get actual code from file
        with open(page., 'r') as f:
            content =

        for line in content:
        if not content or ".. argparse-md::" not in content:
            return None
        # Parse the directive and generate the markdown table
        markdown_table = self.generate_argparse_markdown_table(content)
        if not markdown_table:
            return None

        # Replace the directive with the generated markdown table
        content = content.replace(".. argparse-md::", markdown_table)

        return content

    def generate_argparse_markdown_table(self, content):
        """ Parses the content for the argparse arguments and generates a markdown table"""
        arguments = self.parse_argparse_arguments(content)
        markdown_table = self.generate_markdown_table(arguments)
        return markdown_table

    def parse_argparse_arguments(self, content):
        """ Parses the content for the argparse arguments """
        print(content)
        arguments = []
        for line in content.splitlines():
            cleaned_line = line.strip()
            if cleaned_line.startswith("parser.add_argument"):
                arg_string = cleaned_line[cleaned_line.find(
                    "(")+1:cleaned_line.find(")")]
                arg_list = arg_string.split(",")
                arg_list = [arg.strip() for arg in arg_list]
                for arg in arg_list:
                    arg_dict = dict()
                    for elem in ["help", "type", "default"]:
                        arg_dict[elem] = arg[arg.find(elem):].split("=")[
                            1].strip()
                    arguments.append(arg_dict)
        return arguments

    def generate_markdown_table(self, arguments):
        """ Generates a markdown table from the argparse arguments """
        markdown_table = "| Argument | Type | Default | Help |\n"
        markdown_table += "| -------- | ---- | ------- | ---- |\n"
        for arg in arguments:
            markdown_table += f"| {arg['help']} | {arg['type']} | {arg['default']} | {arg['help']} |\n"
        return markdown_table
