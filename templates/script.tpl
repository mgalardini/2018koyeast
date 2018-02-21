{%- extends 'null.tpl' -%}

{%- block header -%}
#!/usr/bin/env python
# coding: utf-8


def get_options():
    import argparse

    description = ''
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('name',
                        help='Name')

    return parser.parse_args()


if __name__ == "__main__":
    options = get_options()

{% endblock header %}

{% block in_prompt %}
{% endblock in_prompt %}

{% block input %}
{{ cell.source | ipython2python | indent}}
{% endblock input %}

{% block markdowncell scoped %}
{{ cell.source | comment_lines | indent}}
{% endblock markdowncell %}
